from collections import Counter
from pathlib import Path

import parasail


REFERENCE_FASTA = Path(
    "/mnt/share/bhaskar/rcc_data_extend/tally_1k/rcc/"
    "additional_validation/parasail/trna_ref.fa"
)

FASTQ_DIR = Path(
    "/mnt/share/bhaskar/rcc_data_extend/tally_1k/rcc/"
    "fastq_validation"
)

OUTPUT_TSV = Path(
    "/mnt/share/bhaskar/rcc_data_extend/tally_1k/rcc/"
    "additional_validation/parasail/"
    "parasail_best_match_counts.tsv"
)

SUMMARY_OUTPUT_TSV = Path(
    "/mnt/share/bhaskar/rcc_data_extend/tally_1k/rcc/"
    "additional_validation/parasail/"
    "parasail_correct_incorrect_multimap.tsv"
)

# Alignment parameters
GAP_OPEN = 5
GAP_EXTEND = 1
MATRIX = parasail.dnafull


def read_fasta(fasta_path):
    """
    Read reference sequences from a FASTA file.
    """
    references = {}

    current_name = None
    sequence_parts = []

    with open(fasta_path, "r") as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if current_name is not None:
                    references[current_name] = "".join(
                        sequence_parts
                    ).upper()

                current_name = line[1:].split()[0]
                sequence_parts = []

            else:
                sequence_parts.append(line)

    if current_name is not None:
        references[current_name] = "".join(
            sequence_parts
        ).upper()

    return references


def read_fastq_sequences(fastq_path):
    """
    Yield read IDs and sequences from a FASTQ file.
    """
    with open(fastq_path, "r") as handle:
        while True:
            header = handle.readline()

            if not header:
                break

            sequence = handle.readline()
            plus_line = handle.readline()
            quality = handle.readline()

            if not sequence or not plus_line or not quality:
                raise ValueError(
                    f"Incomplete FASTQ record found in {fastq_path}"
                )

            header = header.strip()
            sequence = sequence.strip().upper()

            if not header.startswith("@"):
                raise ValueError(
                    f"Invalid FASTQ header in {fastq_path}: "
                    f"{header}"
                )

            read_id = header[1:].split()[0]

            yield read_id, sequence


def sanitize_sequence(sequence):
    """
    Replace unsupported nucleotide symbols with N.
    """
    allowed_bases = {"A", "C", "G", "T", "N"}

    return "".join(
        base if base in allowed_bases else "N"
        for base in sequence.upper()
    )


def align_read_to_reference(
    read_sequence,
    reference_sequence,
):
    """
    Perform semi-global alignment using parasail.
    """
    result = parasail.sg_striped_16(
        read_sequence,
        reference_sequence,
        GAP_OPEN,
        GAP_EXTEND,
        MATRIX,
    )

    return int(result.score)


def get_best_reference(
    read_sequence,
    references,
):
    """
    Align one read against all references.

    A uniquely highest-scoring reference is returned as a unique match.
    If multiple references share the highest score, the read is
    classified as MULTIMAP.
    """
    best_score = None
    best_references = []

    for reference_name, reference_sequence in references.items():
        score = align_read_to_reference(
            read_sequence=read_sequence,
            reference_sequence=reference_sequence,
        )

        if best_score is None or score > best_score:
            best_score = score
            best_references = [reference_name]

        elif score == best_score:
            best_references.append(reference_name)

    if len(best_references) == 1:
        return (
            best_references[0],
            best_score,
            False,
            1,
        )

    return (
        "MULTIMAP",
        best_score,
        True,
        len(best_references),
    )


def count_best_matches(
    fastq_path,
    references,
):
    """
    Count unique best matches and multimapping reads.
    """
    counts = Counter()

    total_reads = 0
    multimap_reads = 0
    total_multimap_candidates = 0
    total_best_score = 0

    for read_id, read_sequence in read_fastq_sequences(
        fastq_path
    ):
        read_sequence = sanitize_sequence(
            read_sequence
        )

        if not read_sequence:
            counts["EMPTY_SEQUENCE"] += 1
            continue

        (
            best_reference,
            best_score,
            is_multimap,
            number_of_best_references,
        ) = get_best_reference(
            read_sequence=read_sequence,
            references=references,
        )

        counts[best_reference] += 1
        total_reads += 1
        total_best_score += best_score

        if is_multimap:
            multimap_reads += 1
            total_multimap_candidates += (
                number_of_best_references
            )

    return (
        counts,
        total_reads,
        multimap_reads,
        total_multimap_candidates,
        total_best_score,
    )


def normalize_name(name):
    """
    Normalize reference and FASTQ names for comparison.
    """
    return "".join(
        character.lower()
        for character in name
        if character.isalnum()
    )


def get_expected_references(
    fastq_name,
    reference_names,
):
    """
    Identify the expected reference for one FASTQ file.
    """
    expected_name = normalize_name(
        Path(fastq_name).stem
    )

    expected_references = set()

    for reference_name in reference_names:
        if normalize_name(reference_name) == expected_name:
            expected_references.add(reference_name)

    return expected_references


def write_detailed_results(
    all_results,
    reference_names,
    output_path,
):
    """
    Write reference-level assignment counts for each FASTQ file.
    """
    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    columns = (
        [
            "fastq",
            "expected_class",
            "total_reads",
            "unique_reads",
            "multimap_reads",
            "multimap_rate",
            "mean_multimap_candidates",
            "mean_best_score",
        ]
        + reference_names
        + [
            "MULTIMAP",
            "EMPTY_SEQUENCE",
        ]
    )

    with open(output_path, "w") as output_handle:
        output_handle.write(
            "\t".join(columns) + "\n"
        )

        for result in all_results:
            fastq_name = result["fastq_name"]
            counts = result["counts"]
            total_reads = result["total_reads"]
            multimap_reads = result["multimap_reads"]
            total_multimap_candidates = result[
                "total_multimap_candidates"
            ]
            total_best_score = result[
                "total_best_score"
            ]

            expected_class = Path(
                fastq_name
            ).stem

            unique_reads = (
                total_reads
                - multimap_reads
            )

            multimap_rate = (
                multimap_reads / total_reads
                if total_reads > 0
                else 0.0
            )

            mean_multimap_candidates = (
                total_multimap_candidates / multimap_reads
                if multimap_reads > 0
                else 0.0
            )

            mean_best_score = (
                total_best_score / total_reads
                if total_reads > 0
                else 0.0
            )

            row = [
                fastq_name,
                expected_class,
                str(total_reads),
                str(unique_reads),
                str(multimap_reads),
                f"{multimap_rate:.6f}",
                f"{mean_multimap_candidates:.6f}",
                f"{mean_best_score:.6f}",
            ]

            for reference_name in reference_names:
                row.append(
                    str(
                        counts.get(
                            reference_name,
                            0,
                        )
                    )
                )

            row.append(
                str(
                    counts.get(
                        "MULTIMAP",
                        0,
                    )
                )
            )

            row.append(
                str(
                    counts.get(
                        "EMPTY_SEQUENCE",
                        0,
                    )
                )
            )

            output_handle.write(
                "\t".join(row) + "\n"
            )

def write_summary_results(
    all_results,
    reference_names,
    output_path,
):
    """
    Write simple FASTQ-level summary counts.

    Output columns:
    - FASTQ file
    - Best unique reference
    - Best unique count
    - Multimap count
    - Other unique count
    """
    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    columns = [
        "fastq",
        "best_unique_reference",
        "best_unique_count",
        "multimap_count",
        "other_unique_count",
    ]

    with open(output_path, "w") as output_handle:
        output_handle.write(
            "\t".join(columns) + "\n"
        )

        for result in all_results:
            fastq_name = result["fastq_name"]
            counts = result["counts"]

            unique_counts = {
                reference_name: counts.get(
                    reference_name,
                    0,
                )
                for reference_name in reference_names
            }

            if unique_counts:
                best_unique_reference = max(
                    unique_counts,
                    key=unique_counts.get,
                )

                best_unique_count = unique_counts[
                    best_unique_reference
                ]
            else:
                best_unique_reference = ""
                best_unique_count = 0

            total_unique_count = sum(
                unique_counts.values()
            )

            multimap_count = counts.get(
                "MULTIMAP",
                0,
            )

            other_unique_count = (
                total_unique_count
                - best_unique_count
            )

            row = [
                fastq_name,
                best_unique_reference,
                str(best_unique_count),
                str(multimap_count),
                str(other_unique_count),
            ]

            output_handle.write(
                "\t".join(row) + "\n"
            )



def main():
    """
    Run parasail classification for all FASTQ files.
    """
    if not REFERENCE_FASTA.exists():
        raise FileNotFoundError(
            f"Reference FASTA was not found: "
            f"{REFERENCE_FASTA}"
        )

    if not FASTQ_DIR.exists():
        raise FileNotFoundError(
            f"FASTQ directory was not found: "
            f"{FASTQ_DIR}"
        )

    references = read_fasta(
        REFERENCE_FASTA
    )

    if not references:
        raise ValueError(
            f"No reference sequences were found in "
            f"{REFERENCE_FASTA}"
        )

    reference_names = sorted(
        references.keys()
    )

    fastq_paths = sorted(
        FASTQ_DIR.glob("*.fastq")
    )

    if not fastq_paths:
        raise FileNotFoundError(
            f"No FASTQ files were found in "
            f"{FASTQ_DIR}"
        )

    print(
        f"Reference FASTA: {REFERENCE_FASTA}"
    )
    print(
        f"Number of references: {len(reference_names)}"
    )
    print(
        f"FASTQ directory: {FASTQ_DIR}"
    )
    print(
        f"Number of FASTQ files: {len(fastq_paths)}"
    )

    all_results = []

    for fastq_index, fastq_path in enumerate(
        fastq_paths,
        start=1,
    ):
        print()
        print(
            f"[{fastq_index}/{len(fastq_paths)}] "
            f"Processing {fastq_path.name}"
        )

        (
            counts,
            total_reads,
            multimap_reads,
            total_multimap_candidates,
            total_best_score,
        ) = count_best_matches(
            fastq_path=fastq_path,
            references=references,
        )

        all_results.append(
            {
                "fastq_name": fastq_path.name,
                "counts": counts,
                "total_reads": total_reads,
                "multimap_reads": multimap_reads,
                "total_multimap_candidates": (
                    total_multimap_candidates
                ),
                "total_best_score": total_best_score,
            }
        )

        unique_reads = (
            total_reads
            - multimap_reads
        )

        multimap_rate = (
            multimap_reads / total_reads
            if total_reads > 0
            else 0.0
        )

        mean_best_score = (
            total_best_score / total_reads
            if total_reads > 0
            else 0.0
        )

        print(
            f"  Total reads: {total_reads}"
        )
        print(
            f"  Unique best matches: {unique_reads}"
        )
        print(
            f"  Multimap reads: {multimap_reads}"
        )
        print(
            f"  Multimap rate: {multimap_rate:.2%}"
        )
        print(
            f"  Mean best score: {mean_best_score:.2f}"
        )

        print(
            "  Top unique assignments:"
        )

        unique_counts = {
            reference_name: count
            for reference_name, count in counts.items()
            if reference_name not in {
                "MULTIMAP",
                "EMPTY_SEQUENCE",
            }
        }

        for reference_name, count in Counter(
            unique_counts
        ).most_common(5):
            fraction = (
                count / total_reads
                if total_reads > 0
                else 0.0
            )

            print(
                f"    {reference_name}: "
                f"{count} "
                f"({fraction:.2%})"
            )

    write_detailed_results(
        all_results=all_results,
        reference_names=reference_names,
        output_path=OUTPUT_TSV,
    )

    write_summary_results(
        all_results=all_results,
        reference_names=reference_names,
        output_path=SUMMARY_OUTPUT_TSV,
    )

    print()
    print(
        f"Detailed results written to: "
        f"{OUTPUT_TSV}"
    )
    print(
        f"Summary results written to: "
        f"{SUMMARY_OUTPUT_TSV}"
    )


if __name__ == "__main__":
    main()