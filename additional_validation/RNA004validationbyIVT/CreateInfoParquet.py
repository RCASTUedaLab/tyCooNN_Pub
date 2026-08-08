from Bio import SeqIO
from collections import defaultdict
import os
import numpy as np
import pandas as pd
import pysam
import pod5
import pyarrow as pa
import pyarrow.parquet as pq

MAX_PER_TRNA = 12000
MAX_BAM_CANDIDATES_PER_TRNA = 14000

CONSTRUCT_FASTA = "/mnt/share/bhaskar/work/tyCooNN_2206_out/rna004_ivt_test/ivt_batch3.fa"
TRNA_FASTA = "/share/reference/trna/ecolitRNA_unmod.fa"
BAM_PATH = "/mnt/share/ueda/tycoon_revise/ivt_batch3_sup.sorted.bam"
POD5_PATH = "/mnt/share/ueda/tycoon_revise/ivt_batch3.pod5"
OUTPUT_DIR = "/mnt/share/ueda/tycoon_revise/rna004_training"
OUTPUT_REGION_TSV = os.path.join(OUTPUT_DIR, "ivt_batch3_trna_regions.tsv")
OUTPUT_PARQUET = os.path.join(OUTPUT_DIR, "rna004_ivt_training.parquet")
FIXED_LENGTH = 8192

TARGET_TRNAS = [
    "Ser5", "Thr3", "Arg2", "Leu3", "Leu1_P", "Pro3", "Asp", "Val2A",
    "Val2B", "Phe", "Thr2", "Thr4", "Val1", "Gln1", "Ser2", "Arg5",
]


def normalize_sequence(seq):
    return str(seq).upper().replace("U", "T").replace(" ", "").replace("\n", "")


def load_construct(path):
    records = list(SeqIO.parse(path, "fasta"))
    if len(records) != 1:
        raise ValueError("Construct FASTA must contain exactly one sequence.")
    return normalize_sequence(records[0].seq)


def load_trna_sequences(path):
    return {r.id: normalize_sequence(r.seq) for r in SeqIO.parse(path, "fasta")}


def find_trna_sequence(name, seqs):
    if name in seqs:
        return seqs[name]
    hits = [x for x in seqs if name.lower() in x.lower()]
    if len(hits) != 1:
        raise ValueError("Could not uniquely find tRNA: " + name)
    return seqs[hits[0]]


def create_region_table(construct, seqs):
    rows = []
    for trna in TARGET_TRNAS:
        seq = find_trna_sequence(trna, seqs)
        start = construct.find(seq)
        if start < 0:
            raise ValueError("tRNA not found in construct: " + trna)
        rows.append({"trna": trna, "start": start, "end": start + len(seq), "length": len(seq)})
    return pd.DataFrame(rows)


def get_move_info(aln):
    if not aln.has_tag("mv"):
        return None
    mv = list(aln.get_tag("mv"))
    if len(mv) < 2:
        return None
    stride = int(mv[0])
    moves = np.asarray(mv[1:], dtype=np.int16)
    ts = int(aln.get_tag("ts")) if aln.has_tag("ts") else 0
    ns = int(aln.get_tag("ns")) if aln.has_tag("ns") else None
    return stride, moves, ts, ns


def build_base_to_move(moves):
    out = []
    for block, move in enumerate(moves):
        move = int(move)
        if move > 0:
            out.extend([block] * move)
    return np.asarray(out, dtype=np.int64)


def get_query_interval(aln, ref_start, ref_end):
    qpos = []
    for q, r in aln.get_aligned_pairs(matches_only=False):
        if q is None or r is None:
            continue
        if ref_start <= r < ref_end:
            qpos.append(int(q))
    if not qpos:
        return None
    return min(qpos), max(qpos) + 1


def query_to_signal(aln, q_start, q_end, base_to_move, stride, ts):
    qlen = aln.query_length
    if qlen is None:
        return None

    # Direct RNA sequence coordinates are opposite to raw-signal direction.
    b_start = qlen - q_end
    b_end = qlen - q_start

    if b_start < 0 or b_end > len(base_to_move) or b_end <= b_start:
        return None

    move_start = int(base_to_move[b_start])
    move_end = int(base_to_move[b_end - 1] + 1)
    signal_start = int(ts + move_start * stride)
    signal_end = int(ts + move_end * stride)
    return move_start, move_end, signal_start, signal_end


def process_alignment(aln, region_df):
    if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
        return []

    move_info = get_move_info(aln)
    if move_info is None:
        return []

    stride, moves, ts, ns = move_info
    base_to_move = build_base_to_move(moves)
    if len(base_to_move) == 0:
        return []

    rows = []
    for _, region in region_df.iterrows():
        ref_start = int(region["start"])
        ref_end = int(region["end"])

        if aln.reference_end is None:
            continue
        if aln.reference_end <= ref_start or aln.reference_start >= ref_end:
            continue

        qint = get_query_interval(aln, ref_start, ref_end)
        if qint is None:
            continue
        q_start, q_end = qint

        sint = query_to_signal(aln, q_start, q_end, base_to_move, stride, ts)
        if sint is None:
            continue
        move_start, move_end, signal_start, signal_end = sint

        rows.append({
            "read_id": aln.query_name,
            "trna": str(region["trna"]),
            "reference_start": ref_start,
            "reference_end": ref_end,
            "read_start": q_start,
            "read_end": q_end,
            "move_start": move_start,
            "move_end": move_end,
            "stride": stride,
            "ts": ts,
            "ns": ns,
            "signal_start": signal_start,
            "signal_end": signal_end,
            "mapq": int(aln.mapping_quality),
        })
    return rows


def scan_bam(
    path,
    region_df,
):
    rows_by_read = defaultdict(
        list
    )

    n = 0
    n_segments = 0

    trna_counts = {
        trna: 0
        for trna in TARGET_TRNAS
    }

    completed_trnas = set()

    with pysam.AlignmentFile(
        path,
        "rb",
    ) as bam:

        for aln in bam.fetch(
            until_eof=True
        ):

            n += 1

            rows = process_alignment(
                aln,
                region_df,
            )

            for row in rows:

                trna = row[
                    "trna"
                ]

                # Do not keep additional candidates
                # after this class has enough reads.
                if trna in completed_trnas:
                    continue

                read_id = row[
                    "read_id"
                ]

                rows_by_read[
                    read_id
                ].append(
                    row
                )

                trna_counts[
                    trna
                ] += 1

                n_segments += 1

                if (
                    trna_counts[
                        trna
                    ]
                    >= MAX_BAM_CANDIDATES_PER_TRNA
                ):

                    completed_trnas.add(
                        trna
                    )

                    print(
                        "BAM candidate completed:",
                        trna,
                        trna_counts[
                            trna
                        ],
                        "completed:",
                        len(
                            completed_trnas
                        ),
                        "/",
                        len(
                            TARGET_TRNAS
                        ),
                    )

            # Stop BAM scanning as soon as all 16 classes
            # have enough candidate reads.
            if (
                len(
                    completed_trnas
                )
                == len(
                    TARGET_TRNAS
                )
            ):

                print()
                print(
                    "All tRNAs have at least",
                    MAX_BAM_CANDIDATES_PER_TRNA,
                    "BAM candidates."
                )

                print(
                    "Stopping BAM scan at alignment:",
                    n,
                )

                break

            if n % 1000 == 0:

                print()
                print(
                    "BAM alignments:",
                    n,
                    "segments:",
                    n_segments,
                    "reads:",
                    len(
                        rows_by_read
                    ),
                    "completed:",
                    len(
                        completed_trnas
                    ),
                )

                print(
                    "tRNA candidate counts:"
                )

                for trna in TARGET_TRNAS:

                    print(
                        "  %-7s %8d"
                        % (
                            trna,
                            trna_counts[
                                trna
                            ],
                        )
                    )

    print()
    print(
        "Final BAM candidate counts:"
    )

    for trna in TARGET_TRNAS:

        print(
            "  %-7s %8d"
            % (
                trna,
                trna_counts[
                    trna
                ],
            )
        )

    return rows_by_read


def median_mad_normalize(signal):
    signal = np.asarray(signal, dtype=np.float32)
    med = float(np.median(signal))
    mad = float(np.median(np.abs(signal - med)))
    scale = 1.4826 * mad
    if scale < 1e-6:
        scale = float(np.std(signal))
    if scale < 1e-6:
        scale = 1.0
    return (signal - med) / scale


def pad_signal(signal, fixed_length):
    out = np.zeros(fixed_length, dtype=np.float32)
    n = min(len(signal), fixed_length)
    out[:n] = signal[:n]
    return out

def write_training_parquet(
    pod5_path,
    rows_by_read,
    output_path,
    fixed_length,
):
    schema = pa.schema([
        ("read_id", pa.string()),
        ("trna", pa.string()),
        ("trimsignal", pa.list_(pa.float32())),
        ("raw_signal_length", pa.int32()),
        ("reference_start", pa.int32()),
        ("reference_end", pa.int32()),
        ("read_start", pa.int32()),
        ("read_end", pa.int32()),
        ("move_start", pa.int32()),
        ("move_end", pa.int32()),
        ("signal_start", pa.int64()),
        ("signal_end", pa.int64()),
        ("mapq", pa.int32()),
    ])

    writer = pq.ParquetWriter(
        output_path,
        schema,
        compression="snappy",
    )

    buffer = []

    n_reads = 0
    n_segments = 0
    n_rejected = 0

    # Count accepted signals for each tRNA.
    trna_counts = {
        trna: 0
        for trna in TARGET_TRNAS
    }

    # Keep track of tRNAs that reached MAX_PER_TRNA.
    completed_trnas = set()

    with pod5.Reader(
        pod5_path
    ) as reader:

        for read in reader.reads():

            # Stop POD5 scanning immediately when all classes are complete.
            if len(completed_trnas) == len(TARGET_TRNAS):
                print(
                    "All tRNAs completed:",
                    len(completed_trnas),
                )
                break

            read_id = str(
                read.read_id
            )

            if read_id not in rows_by_read:
                continue

            n_reads += 1

            signal_pa = np.asarray(
                read.signal_pa,
                dtype=np.float32,
            )

            for row in rows_by_read[
                read_id
            ]:

                trna = row[
                    "trna"
                ]

                # Skip classes that are already complete.
                if trna in completed_trnas:
                    continue

                s = int(
                    row[
                        "signal_start"
                    ]
                )

                e = int(
                    row[
                        "signal_end"
                    ]
                )

                if (
                    s < 0
                    or e > len(signal_pa)
                    or e <= s
                ):
                    n_rejected += 1
                    continue

                raw = signal_pa[
                    s:e
                ]

                if len(raw) == 0:
                    n_rejected += 1
                    continue

                norm = (
                    median_mad_normalize(
                        raw
                    )
                )

                fixed = pad_signal(
                    norm,
                    fixed_length,
                )

                buffer.append(
                    {
                        "read_id":
                            read_id,

                        "trna":
                            trna,

                        "trimsignal":
                            fixed.tolist(),

                        "raw_signal_length":
                            int(
                                len(raw)
                            ),

                        "reference_start":
                            int(
                                row[
                                    "reference_start"
                                ]
                            ),

                        "reference_end":
                            int(
                                row[
                                    "reference_end"
                                ]
                            ),

                        "read_start":
                            int(
                                row[
                                    "read_start"
                                ]
                            ),

                        "read_end":
                            int(
                                row[
                                    "read_end"
                                ]
                            ),

                        "move_start":
                            int(
                                row[
                                    "move_start"
                                ]
                            ),

                        "move_end":
                            int(
                                row[
                                    "move_end"
                                ]
                            ),

                        "signal_start":
                            s,

                        "signal_end":
                            e,

                        "mapq":
                            int(
                                row[
                                    "mapq"
                                ]
                            ),
                    }
                )

                trna_counts[
                    trna
                ] += 1

                n_segments += 1

                # Mark this tRNA as completed immediately.
                if (
                    trna_counts[
                        trna
                    ]
                    >= MAX_PER_TRNA
                ):

                    completed_trnas.add(
                        trna
                    )

                    print(
                        "Completed:",
                        trna,
                        trna_counts[
                            trna
                        ],
                        "completed classes:",
                        len(
                            completed_trnas
                        ),
                        "/",
                        len(
                            TARGET_TRNAS
                        ),
                    )

                if len(buffer) >= 1000:

                    writer.write_table(
                        pa.Table.from_pylist(
                            buffer,
                            schema=schema,
                        )
                    )

                    buffer = []

            # Break immediately after finishing the current read
            # if all 16 classes have reached MAX_PER_TRNA.
            if len(completed_trnas) == len(TARGET_TRNAS):

                print(
                    "Collected",
                    MAX_PER_TRNA,
                    "signals for all",
                    len(TARGET_TRNAS),
                    "tRNAs."
                )

                break

            if (
                n_reads
                % 1000
                == 0
            ):

                print(
                    "POD5 reads:",
                    n_reads,
                    "written:",
                    n_segments,
                    "rejected:",
                    n_rejected,
                    "completed:",
                    len(
                        completed_trnas
                    ),
                )

                print(
                    trna_counts
                )

    if buffer:

        writer.write_table(
            pa.Table.from_pylist(
                buffer,
                schema=schema,
            )
        )

    writer.close()

    print()
    print(
        "Saved:",
        output_path
    )

    print(
        "Written segments:",
        n_segments
    )

    print(
        "Rejected segments:",
        n_rejected
    )

    print()
    print(
        "Final counts:"
    )

    for trna in TARGET_TRNAS:

        print(
            trna,
            trna_counts[
                trna
            ],
        )


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    construct = load_construct(CONSTRUCT_FASTA)
    trna_sequences = load_trna_sequences(TRNA_FASTA)
    region_df = create_region_table(construct, trna_sequences)

    print(region_df.to_string(index=False))
    region_df.to_csv(OUTPUT_REGION_TSV, sep="\t", index=False)

    rows_by_read = scan_bam(BAM_PATH, region_df)
    write_training_parquet(POD5_PATH, rows_by_read, OUTPUT_PARQUET, FIXED_LENGTH)


if __name__ == "__main__":
    main()
