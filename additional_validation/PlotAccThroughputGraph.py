import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


# ==========================================================
# SVG text setting
# ==========================================================

plt.rcParams["svg.fonttype"] = "none"


# ==========================================================
# Input / output
# ==========================================================

INPUT_TSV = "/mnt/share/bhaskar/rcc_data_extend/tally_1k/rcc/additional_validation/comparison_methods.tsv"
OUTPUT_PNG = "/mnt/share/bhaskar/rcc_data_extend/tally_1k/rcc/additional_validation/comparison_with_parasail_CD.png"
OUTPUT_PDF = "/mnt/share/bhaskar/rcc_data_extend/tally_1k/rcc/additional_validation/comparison_with_parasail_CD.pdf"
OUTPUT_SVG = "/mnt/share/bhaskar/rcc_data_extend/tally_1k/rcc/additional_validation/comparison_with_parasail_CD.svg"


# ==========================================================
# Plot style
# ==========================================================

sns.set_style("white")

METHOD_ORDER = [
    "ntc44",
    "BWA-MEM",
    "Minimap2",
    "parasail",
]

PALETTE = {
    "ntc44": "#E64B35",
    "BWA-MEM": "#4C78A8",
    "Minimap2": "#C9AD3A",
    "parasail": "#7F7F7F",   # green,
}

MARKERS = {
    "ntc44": "o",
    "BWA-MEM": "^",
    "Minimap2": "s",
    "parasail": "D",
}


# ==========================================================
# Utilities
# ==========================================================

def clean_percent(value):
    """
    Convert percentage strings such as '75.5%' to float 75.5.
    """
    if pd.isna(value):
        return np.nan

    if isinstance(value, str):
        value = value.strip().replace("%", "")
        if value == "":
            return np.nan

    return float(value)


def simplify_method_name(method_name):
    """
    Convert long header names into simplified method names.
    """
    method_lower = str(method_name).lower()

    if "minimap2" in method_lower:
        return "Minimap2"
    elif "bwa-mem" in method_lower:
        return "BWA-MEM"
    elif "parasail" in method_lower:
        return "parasail"
    elif "ntc44" in method_lower:
        return "ntc44"
    else:
        return None


def find_metric_column(columns, method_name, keyword):
    """
    Find a sub-column under one method whose second-level header
    contains the given keyword.
    """
    for col in columns:
        if col[0] == method_name and keyword.lower() in str(col[1]).lower():
            return col
    return None


# ==========================================================
# Data loading
# ==========================================================

def load_method_comparison_table(input_tsv):
    """
    Read the two-header TSV table and return a tidy DataFrame.

    Output columns:
      - tRNA_species
      - Method
      - Accuracy
      - Read_retention_rate
    """
    df_raw = pd.read_csv(
        input_tsv,
        sep="\t",
        header=[0, 1],
    )

    cleaned_columns = []
    for col0, col1 in df_raw.columns:
        col0 = str(col0).strip()
        col1 = str(col1).strip()
        cleaned_columns.append((col0, col1))
    df_raw.columns = pd.MultiIndex.from_tuples(cleaned_columns)

    species_col = df_raw.columns[0]
    species_series = df_raw[species_col].astype(str).str.strip()

    keep_mask = (
        species_series.notna()
        & (species_series != "")
        & (~species_series.str.lower().str.startswith("average"))
    )
    df_raw = df_raw.loc[keep_mask].copy()

    tidy_rows = []

    unique_method_headers = list(dict.fromkeys(df_raw.columns.get_level_values(0)))

    for raw_method in unique_method_headers:
        simplified_method = simplify_method_name(raw_method)

        if simplified_method is None:
            continue

        accuracy_col = find_metric_column(
            df_raw.columns,
            raw_method,
            "Accuracy",
        )
        retention_col = find_metric_column(
            df_raw.columns,
            raw_method,
            "Read retention rate",
        )

        if accuracy_col is None or retention_col is None:
            print(f"Warning: Could not find required columns for method {raw_method}")
            continue

        for idx in df_raw.index:
            tidy_rows.append(
                {
                    "tRNA_species": str(df_raw.loc[idx, species_col]).strip(),
                    "Method": simplified_method,
                    "Accuracy": clean_percent(df_raw.loc[idx, accuracy_col]),
                    "Read_retention_rate": clean_percent(df_raw.loc[idx, retention_col]),
                }
            )

    df_tidy = pd.DataFrame(tidy_rows)

    df_tidy["Method"] = pd.Categorical(
        df_tidy["Method"],
        categories=METHOD_ORDER,
        ordered=True,
    )

    return df_tidy


# ==========================================================
# Panel C
# ==========================================================

def draw_panel_c(fig, outer_spec, df_tidy):
    left_spec = outer_spec.subgridspec(
        2,
        2,
        width_ratios=[4.0, 1.0],
        height_ratios=[1.0, 4.0],
        wspace=0.05,
        hspace=0.05,
    )

    ax_top = fig.add_subplot(left_spec[0, 0])
    ax_scatter = fig.add_subplot(left_spec[1, 0])
    ax_right = fig.add_subplot(left_spec[1, 1], sharey=ax_scatter)
    ax_empty = fig.add_subplot(left_spec[0, 1])
    ax_empty.axis("off")

    for method in METHOD_ORDER:
        sub = df_tidy[df_tidy["Method"] == method].copy()

        if sub.empty:
            continue

        ax_scatter.scatter(
            sub["Accuracy"],
            sub["Read_retention_rate"],
            s=14,
            alpha=0.85,
            color=PALETTE[method],
            marker=MARKERS[method],
            label=method,
            edgecolor="none",
        )

        sns.kdeplot(
            data=sub,
            x="Accuracy",
            ax=ax_top,
            color=PALETTE[method],
            lw=0.8,
            fill=False,
            clip=(0, 100),
        )

        sns.kdeplot(
            data=sub,
            y="Read_retention_rate",
            ax=ax_right,
            color=PALETTE[method],
            lw=0.8,
            fill=False,
            clip=(0, 100),
        )

    ax_scatter.set_xlim(0, 100)
    ax_scatter.set_ylim(0, 100)
    ax_scatter.set_xlabel("Accuracy (%)")
    ax_scatter.set_ylabel("Read retention rate (%)")
    ax_scatter.legend(
        frameon=False,
        loc="lower left",
        fontsize=9,
    )

    ax_top.set_xlim(0, 100)
    ax_top.set_xlabel("")
    ax_top.set_ylabel("Density")
    ax_top.tick_params(axis="x", labelbottom=False)

    ax_right.set_ylim(0, 100)
    ax_right.set_xlabel("Density")
    ax_right.set_ylabel("")
    ax_right.tick_params(axis="y", labelleft=False)

    ax_top.text(
        -0.12,
        1.25,
        "c",
        transform=ax_top.transAxes,
        fontsize=18,
        fontweight="bold",
        va="top",
        ha="left",
    )


# ==========================================================
# Panel D
# ==========================================================

def draw_panel_d(fig, outer_spec, df_tidy):
    right_spec = outer_spec.subgridspec(
        2,
        1,
        hspace=0.18,
    )

    ax_acc = fig.add_subplot(right_spec[0, 0])
    ax_ret = fig.add_subplot(right_spec[1, 0], sharex=ax_acc)

    sns.boxplot(
        data=df_tidy,
        x="Method",
        y="Accuracy",
        order=METHOD_ORDER,
        palette=PALETTE,
        linewidth=1.2,
        fliersize=3,
        width=0.7,
        ax=ax_acc,
    )

    ax_acc.set_ylabel("Accuracy (%)")
    ax_acc.set_xlabel("")
    ax_acc.set_ylim(0, 105)
    ax_acc.tick_params(axis="x", bottom=False, labelbottom=False, top=True, labeltop=True)
    ax_acc.set_xticklabels(
        METHOD_ORDER,
        rotation=45,
        ha="left",
    )

    acc_means = (
        df_tidy.groupby("Method", observed=False)["Accuracy"]
        .mean()
        .reindex(METHOD_ORDER)
    )

    for i, method in enumerate(METHOD_ORDER):
        value = acc_means.loc[method]
        if pd.notna(value):
            ax_acc.text(
                i,
                min(value + 4, 102),
                f"{value:.1f}",
                ha="center",
                va="bottom",
                fontsize=10,
            )

    sns.boxplot(
        data=df_tidy,
        x="Method",
        y="Read_retention_rate",
        order=METHOD_ORDER,
        palette=PALETTE,
        linewidth=1.2,
        fliersize=3,
        width=0.7,
        ax=ax_ret,
    )

    ax_ret.set_ylabel("Read retention rate (%)")
    ax_ret.set_xlabel("")
    ax_ret.set_ylim(0, 100)
    ax_ret.set_xticklabels(
        METHOD_ORDER,
        rotation=45,
        ha="right",
    )

    ret_means = (
        df_tidy.groupby("Method", observed=False)["Read_retention_rate"]
        .mean()
        .reindex(METHOD_ORDER)
    )

    for i, method in enumerate(METHOD_ORDER):
        value = ret_means.loc[method]
        if pd.notna(value):
            ax_ret.text(
                i,
                min(value + 4, 98),
                f"{value:.1f}",
                ha="center",
                va="bottom",
                fontsize=10,
            )

    ax_acc.text(
        -0.22,
        1.12,
        "d",
        transform=ax_acc.transAxes,
        fontsize=18,
        fontweight="bold",
        va="top",
        ha="left",
    )


# ==========================================================
# Main
# ==========================================================

def main():
    df_tidy = load_method_comparison_table(INPUT_TSV)

    print(df_tidy.head())
    print()
    print("Mean values by method:")
    print(
        df_tidy.groupby("Method", observed=False)[
            ["Accuracy", "Read_retention_rate"]
        ].mean()
    )
    fig = plt.figure(figsize=(9, 6))

    outer = fig.add_gridspec(
        1,
        2,
        width_ratios=[1.3, 0.6],
        wspace=0.28,
    )


    draw_panel_c(fig, outer[0], df_tidy)
    draw_panel_d(fig, outer[1], df_tidy)

    plt.tight_layout()

    fig.savefig(
        OUTPUT_PNG,
        dpi=300,
        bbox_inches="tight",
    )
    fig.savefig(
        OUTPUT_PDF,
        bbox_inches="tight",
    )
    fig.savefig(
        OUTPUT_SVG,
        bbox_inches="tight",
    )

    plt.show()


if __name__ == "__main__":
    main()