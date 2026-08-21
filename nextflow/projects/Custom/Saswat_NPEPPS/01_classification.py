# 01_classification.py
# ============
# Batch TME (Tumor Microenvironment) classification across multiple cohorts
# using the BostonGene MFP method (Bagaev et al. 2021).

# Each cohort is processed independently:
  # TPM tsv  →  log2(1+x)  →  dedup symbols  →  ssGSEA  →  median_scale  →  KNN predict
  # →  per-cohort TSV  +  per-cohort QC plots

# All 40 per-cohort results are merged into one combined TSV at the end.

# Usage
# -----
    # conda activate omics
    # python 01_classification.py

    # # Skip QC plot generation (faster first-pass):
    # python 01_classification.py --no-qc

# Outputs (written to RESULTS_DIR)
# ---------
    # classified_samples/   one TSV per cohort
    # qc_plots/             one sub-folder per cohort with 4 QC plots
    # combined_classification.tsv   all cohorts stacked, columns: sample_id, TME, cohort
    # run_batch.log         timing + per-dataset stats

# Author : generated for Saswat_NPEPPS project

# ── paths ─────────────────────────────────────────────────────────────────────
# Change these two lines when moving to a new machine/cluster.
PORTRAITS_DIR = "/home/kailasamms/resources/BostonGeneMFP/portraits"
REFERENCE_DIR = "/home/kailasamms/resources/BostonGeneMFP/reference"

#METADATA_DIR = None  # Use None if you want to use all samples
METADATA_DIR = "/home/kailasamms/scripts/nextflow/resources/Datasets/01.Metadata/"
COUNTS_DIR   = "/home/kailasamms/scripts/nextflow/resources/Datasets/02.Counts/"
RESULTS_DIR  = "/home/kailasamms/scripts/nextflow/projects/Custom/Saswat_NPEPPS/01.Classification/"

# ── input file naming convention ────────────────────────────────────────────
# Files must look like: {COHORT}{TPM_SUFFIX}, e.g. TCGA-CHOL_TPM_Counts.tsv
TPM_SUFFIX = "_TPM_Counts.tsv"
METADATA_SUFFIX = "_metadata.xlsx"

# ── imports ───────────────────────────────────────────────────────────────────
import argparse
import logging
import os
import sys
import time
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")          # non-interactive backend — safe on HPC (no display needed)
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Add portraits package (not pip-installable; lives in shared resources folder)
sys.path.insert(0, str(Path(PORTRAITS_DIR).parent))
from portraits.utils import read_gene_sets, ssgsea_formula, median_scale
from portraits.classification import KNeighborsClusterClassifier

warnings.filterwarnings("ignore", category=FutureWarning)


# ── logging ───────────────────────────────────────────────────────────────────
def setup_logging(results_dir: Path) -> logging.Logger:
    results_dir.mkdir(parents=True, exist_ok=True)
    log_path = results_dir / "run_batch.log"

    fmt = "%(asctime)s  %(levelname)-8s  %(message)s"
    logging.basicConfig(
        level=logging.INFO,
        format=fmt,
        handlers=[
            logging.FileHandler(log_path),
            logging.StreamHandler(sys.stdout),
        ],
    )
    return logging.getLogger("run_batch")


# ── cohort naming (once, outside the loop) ────────────────────────────────
def cohort_name_from_path(tsv_path: Path) -> str:
    """Derive cohort name by stripping the configured TPM_SUFFIX from the filename."""
    name = tsv_path.name
    if name.endswith(TPM_SUFFIX):
        return name[: -len(TPM_SUFFIX)]
    # Fallback: suffix didn't match as expected — warn loudly rather than
    # silently producing a mangled name.
    return tsv_path.stem


# ── reference loading (once, outside the loop) ────────────────────────────────
def load_reference(reference_dir: Path, log: logging.Logger):
    """Load TCGA reference signatures, annotations, and fit the KNN model."""
    log.info("Loading TCGA reference cohort from %s", reference_dir)

    tcga_scores = pd.read_csv(reference_dir / "signatures.tsv", sep="\t", index_col=0)
    tcga_annot  = pd.read_csv(reference_dir / "annotation.tsv",  sep="\t", index_col=0)
    gene_sigs   = read_gene_sets(str(reference_dir / "gene_signatures.gmt"))

    model = KNeighborsClusterClassifier(norm=False, scale=False, clip=2, k=35).fit(
        tcga_scores.T, tcga_annot.MFP
    )

    log.info(
        "Reference loaded: %d TCGA samples, %d gene signatures, %d model features",
        tcga_scores.shape[1], len(gene_sigs), len(model.X.columns),
    )
    return model, gene_sigs


# ── per-dataset loader ────────────────────────────────────────────────────────
def load_tpm(tsv_path: Path, log: logging.Logger) -> pd.DataFrame:
    """
    Read a {COHORT}_TPM.tsv file and return a genes-x-samples DataFrame
    of log2(1+TPM) values with gene symbols as the index.

    Steps
    -----
    1. Read the tab-separated file.
    2. Set SYMBOL as index; drop the Ensembl ID and gene_type columns.
    3. Collapse duplicate gene symbols: keep the row with highest total
       expression across all samples (mirrors plot_heatmap() convention).
    4. Log2-transform: always applied to raw TPM, regardless of value range.
    """
    t0 = time.time()

    df = pd.read_csv(tsv_path, sep="\t")
    log.info("  Read %d rows x %d cols in %.1fs", *df.shape, time.time() - t0)

    # ── set SYMBOL as index ───────────────────────────────────────────────────
    if "SYMBOL" not in df.columns:
        raise ValueError(f"No 'SYMBOL' column found in {tsv_path.name}")
    df = df.drop(columns=["ID", "gene_type"], errors="ignore")
    df = df.set_index("SYMBOL")

    n_before = len(df)

    # ── collapse duplicate gene symbols ──────────────────────────────────────
    if df.index.duplicated().any():
        n_dups = df.index.duplicated(keep=False).sum()
        df = df.assign(_total_expr=df.sum(axis=1))
        df = df.sort_values("_total_expr", ascending=False)
        df = df.loc[~df.index.duplicated(keep="first")]   # ← now uses sorted df.index ✓
        df = df.drop(columns="_total_expr")
        n_after = len(df)
        log.info(
            "  Duplicate SYMBOL collapse: %d rows → %d rows "
            "(%d rows involved in duplications, kept highest-expressing copy)",
            n_before, n_after, n_dups,
        )

    # ── log2(1 + TPM) ─────────────────────────────────────────────────────────
    df = np.log2(df + 1)

    log.info(
        "  Final expression matrix: %d genes × %d samples  "
        "value range [%.2f, %.2f]",
        *df.shape, df.min().min(), df.max().max(),
    )
    return df


def filter_tumor_samples_by_column(df: pd.DataFrame, metadata_dir: Path | None,
    cohort: str, log: logging.Logger) -> pd.DataFrame:

    # Metadata filtering is optional
    if metadata_dir is None:
        log.info("  [%s] Metadata filtering disabled. Using all samples.", cohort)
        return df

    metadata_path = metadata_dir / f"{cohort}{METADATA_SUFFIX}"

    if not metadata_path.exists():
        log.info(
            "  [%s] Metadata file not found at %s. Using all samples.",
            cohort, metadata_path
        )
        return df

    #meta_df = pd.read_csv(metadata_path, sep="\t")
    meta_df = pd.read_excel(metadata_path)

    sample_col = "Sample_ID" if "Sample_ID" in meta_df.columns else meta_df.columns[0]

    if "Tissue_Type" not in meta_df.columns:
        log.warning(
            "  [%s] 'Tissue_Type' column not found in metadata. Using all samples.",
            cohort
        )
        return df

    tumor_samples = (
        meta_df.loc[
            meta_df["Tissue_Type"].astype(str).str.strip().str.lower() == "tumor",
            sample_col
        ]
        .astype(str)
        .tolist()
    )

    valid_tumor_cols = [c for c in tumor_samples if c in df.columns]

    filtered_df = df[valid_tumor_cols]

    log.info(
        "  [%s] Tumor filtering: %d/%d tumor samples retained.",
        cohort, len(valid_tumor_cols), len(tumor_samples)
    )

    return filtered_df


# ── gene-coverage check ───────────────────────────────────────────────────────
def check_gene_coverage(expr: pd.DataFrame, gene_sigs: dict,
                        model_cols, log: logging.Logger, cohort: str):
    """
    Warn if a dataset is missing signature genes or model features.
    Neither issue hard-fails; both silently produce zeros/NaN in downstream
    tools, so explicit logging here is the only way to catch it.
    """
    all_sig_genes = set()
    for gs in gene_sigs.values():
        all_sig_genes.update(gs.genes)

    missing_sig = all_sig_genes - set(expr.index)
    pct_missing = 100 * len(missing_sig) / len(all_sig_genes)
    if pct_missing > 0:
        log.warning(
            "[%s] %d / %d signature genes missing (%.1f%%) — "
            "ssGSEA scores will be 0 for affected gene sets",
            cohort, len(missing_sig), len(all_sig_genes), pct_missing,
        )


# ── QC plots ──────────────────────────────────────────────────────────────────
def make_qc_plots(expr: pd.DataFrame, sig_scores_scaled: pd.DataFrame,
                  out_dir: Path, cohort: str, log: logging.Logger):
    """
    Save four QC plots mirroring the original notebook:
        1. Gene expression distribution (all samples overlaid)
        2. PCA — all samples (colour = PCA outlier flag)
        3. Signature score distribution after median scaling
        4. Signature scores PCA
    UMAP is skipped for single-cohort runs (no batch labels available).
    """
    from sklearn.decomposition import PCA
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Expression distribution
    fig, ax = plt.subplots(figsize=(8, 4))
    for col in expr.columns:
        ax.hist(expr[col].dropna(), bins=80, alpha=0.05, color="steelblue", density=True)
    ax.set_xlabel("log2(1+TPM)")
    ax.set_ylabel("Density")
    ax.set_title(f"{cohort} — Expression distribution (all samples)")
    fig.tight_layout()
    fig.savefig(out_dir / "01_expression_distribution.png", dpi=120)
    plt.close(fig)

    # 2. PCA on expression — flag outliers by Mahalanobis-like PCA distance
    try:
        pca = PCA(n_components=2, random_state=0)
        coords = pca.fit_transform(expr.T.dropna(axis=1))
        center = coords.mean(axis=0)
        dist   = np.linalg.norm(coords - center, axis=1)
        threshold = dist.mean() + 3 * dist.std()
        is_outlier = dist > threshold

        fig, ax = plt.subplots(figsize=(7, 6))
        ax.scatter(coords[~is_outlier, 0], coords[~is_outlier, 1],
                   s=15, alpha=0.6, label="Normal")
        if is_outlier.any():
            ax.scatter(coords[is_outlier, 0], coords[is_outlier, 1],
                       s=30, alpha=0.9, color="red", label=f"Outlier (n={is_outlier.sum()})")
        ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
        ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
        ax.set_title(f"{cohort} — PCA (expression)")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(out_dir / "02_pca_expression.png", dpi=120)
        plt.close(fig)

        n_outliers = int(is_outlier.sum())
        if n_outliers > 0:
            log.warning("[%s] PCA flagged %d potential outlier sample(s) — see qc_plots",
                        cohort, n_outliers)
    except Exception as e:
        log.warning("[%s] PCA plot failed: %s", cohort, e)

    # 3. Signature score distribution after scaling
    fig, ax = plt.subplots(figsize=(8, 4))
    for col in sig_scores_scaled.columns:
        ax.hist(sig_scores_scaled[col].dropna(), bins=40,
                alpha=0.4, density=True, label=col)
    ax.set_xlabel("Median-scaled ssGSEA score")
    ax.set_ylabel("Density")
    ax.set_title(f"{cohort} — Scaled signature score distribution")
    ax.legend(fontsize=6, ncol=3)
    fig.tight_layout()
    fig.savefig(out_dir / "03_signature_score_distribution.png", dpi=120)
    plt.close(fig)

    # 4. PCA on signature scores
    try:
        valid = sig_scores_scaled.dropna()
        if len(valid) >= 3:
            pca2 = PCA(n_components=2, random_state=0)
            coords2 = pca2.fit_transform(valid)
            fig, ax = plt.subplots(figsize=(7, 6))
            ax.scatter(coords2[:, 0], coords2[:, 1], s=15, alpha=0.6)
            ax.set_xlabel(f"PC1 ({pca2.explained_variance_ratio_[0]*100:.1f}%)")
            ax.set_ylabel(f"PC2 ({pca2.explained_variance_ratio_[1]*100:.1f}%)")
            ax.set_title(f"{cohort} — PCA (signature scores)")
            fig.tight_layout()
            fig.savefig(out_dir / "04_pca_signature_scores.png", dpi=120)
            plt.close(fig)
    except Exception as e:
        log.warning("[%s] Signature PCA plot failed: %s", cohort, e)

    log.info("[%s] QC plots saved to %s", cohort, out_dir)


# ── classify one dataset ──────────────────────────────────────────────────────
def classify_cohort(tsv_path: Path, model, gene_sigs,
                    results_dir: Path, run_qc: bool,
                    log: logging.Logger) -> pd.Series | None:
    """
    Full pipeline for one cohort.
    Returns a pd.Series (index=sample_id, values=TME) or None on failure.
    """
    cohort = cohort_name_from_path(tsv_path)
    log.info("=" * 60)
    log.info("Processing: %s", cohort)
    t_start = time.time()

    try:
        # 1. Load + normalize
        t0 = time.time()
        expr = load_tpm(tsv_path, log)
        log.info("  [%s] Load+normalize: %.1fs", cohort, time.time() - t0)

        # Filter if metadata provided
        expr = filter_tumor_samples_by_column(expr, Path(METADATA_DIR) if METADATA_DIR else None, cohort, log)

        # 2. Coverage check
        check_gene_coverage(expr, gene_sigs, model.X.columns, log, cohort)

        # 3. ssGSEA
        t0 = time.time()
        # ssgsea_formula expects samples in rows, genes in columns
        sig_scores = ssgsea_formula(expr.T, gene_sigs)
        log.info("  [%s] ssGSEA: %.1fs", cohort, time.time() - t0)

        # 4. Median scale
        sig_scores_scaled = median_scale(sig_scores)

        # 4b. Save ssGSEA scores  ← ADD THIS
        ssgsea_dir = results_dir / "ssgsea_scores"
        ssgsea_dir.mkdir(parents=True, exist_ok=True)
        ssgsea_path = ssgsea_dir / f"{cohort}_ssgsea_scores.tsv"
        sig_scores_scaled.index.name = "sample_id"
        sig_scores_scaled.to_csv(ssgsea_path, sep="\t")
        log.info("  [%s] ssGSEA scores saved → %s", cohort, ssgsea_path)
        
        # 5. Classify
        t0 = time.time()
        valid_scores = sig_scores_scaled[model.X.columns].dropna()
        n_dropped = len(sig_scores_scaled) - len(valid_scores)
        if n_dropped > 0:
            log.warning("[%s] %d sample(s) dropped before predict (NaN in model features)",
                        cohort, n_dropped)

        labels = model.predict(valid_scores).rename("TME")
        log.info("  [%s] Classification: %.1fs  →  %d samples classified",
                 cohort, time.time() - t0, len(labels))
        log.info("  [%s] TME breakdown: %s",
                 cohort, labels.value_counts().to_dict())

        # 6. Save per-cohort TSV
        classified_dir = results_dir / "classified_samples"
        classified_dir.mkdir(parents=True, exist_ok=True)
        out_path = classified_dir / f"{cohort}_classified.tsv"
        labels.to_csv(out_path, sep="\t", header=True)
        log.info("  [%s] Saved → %s", cohort, out_path)

        # 7. QC plots
        if run_qc:
            t0 = time.time()
            qc_dir = results_dir / "qc_plots" / cohort
            make_qc_plots(expr, sig_scores_scaled, qc_dir, cohort, log)
            log.info("  [%s] QC plots: %.1fs", cohort, time.time() - t0)

        log.info("[%s] Total: %.1fs", cohort, time.time() - t_start)
        return labels.rename(cohort)     # carry cohort name for combining

    except Exception as e:
        log.error("[%s] FAILED — %s", cohort, e, exc_info=True)
        return None


# ── main ──────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="Batch TME classification")
    parser.add_argument("--no-qc", action="store_true",
                        help="Skip QC plot generation (faster)")
    args = parser.parse_args()
    run_qc = not args.no_qc

    results_dir = Path(RESULTS_DIR)
    log = setup_logging(results_dir)

    log.info("run_batch.py starting")
    log.info("PORTRAITS_DIR : %s", PORTRAITS_DIR)
    log.info("REFERENCE_DIR : %s", REFERENCE_DIR)
    log.info("COUNTS_DIR    : %s", COUNTS_DIR)
    log.info("RESULTS_DIR   : %s", RESULTS_DIR)
    log.info("QC plots      : %s", "enabled" if run_qc else "disabled (--no-qc)")

    # ── load reference once ───────────────────────────────────────────────────
    model, gene_sigs = load_reference(Path(REFERENCE_DIR), log)

    # ── discover input files ──────────────────────────────────────────────────
    counts_dir = Path(COUNTS_DIR)
    tsv_files = sorted(counts_dir.glob(f"*{TPM_SUFFIX}"))
    if not tsv_files:
        log.error("No *{TPM_SUFFIX} files found in %s — nothing to do.", counts_dir)
        sys.exit(1)
    log.info("Found %d cohort(s): %s",
             len(tsv_files), [f.stem for f in tsv_files])

    # ── save BostonGene 22-pathway gene set reference xlsx ────────────────────
    HEATMAP_22 = [
        "MHCI", "MHCII", "Coactivation_molecules", "Effector_cells",
        "T_cell_traffic", "NK_cells", "T_cells", "B_cells",
        "M1_signatures", "Th1_signature", "Antitumor_cytokines",
        "Checkpoint_inhibition", "Treg", "T_reg_traffic",
        "Neutrophil_signature", "Granulocyte_traffic", "MDSC",
        "MDSC_traffic", "Macrophages", "Macrophage_DC_traffic",
        "Th2_signature", "Protumor_cytokines",
    ]

    try:
        import openpyxl
        from openpyxl import Workbook
        wb = Workbook()
        ws = wb.active
        ws.title = "BostonGene_22_Pathways"
        ws.append(["Pathway"] + [f"Gene_{i+1}" for i in range(20)])  # header
        for sig_name in HEATMAP_22:
            sig = gene_sigs[sig_name]
            genes = list(sig.genes) if hasattr(sig, "genes") else list(sig)
            ws.append([sig_name] + genes)
        ref_path = results_dir / "BostonGene_22_pathway_gene_sets.xlsx"
        wb.save(ref_path)
        log.info("BostonGene 22-pathway reference saved → %s", ref_path)
    except Exception as e:
        log.warning("Could not save pathway reference xlsx: %s", e)

    # ── loop over cohorts ─────────────────────────────────────────────────────
    all_labels = []
    failed = []

    for tsv_path in tsv_files:
        cohort = cohort_name_from_path(tsv_path)
        labels = classify_cohort(
            tsv_path, model, gene_sigs, results_dir, run_qc, log
        )
        if labels is not None:
            all_labels.append(
                pd.DataFrame({
                    "sample_id": labels.index,
                    "TME":       labels.values,
                    "cohort":    cohort,
                })
            )
        else:
            failed.append(cohort)

    # ── combine ssGSEA scores ─────────────────────────────────────────────────
    ssgsea_parts = []
    for tsv_path in tsv_files:
        cohort = cohort_name_from_path(tsv_path)
        score_file = results_dir / "ssgsea_scores" / f"{cohort}_ssgsea_scores.tsv"
        if score_file.exists():
            df_s = pd.read_csv(score_file, sep="\t", index_col=0)
            df_s.insert(0, "cohort", cohort)
            ssgsea_parts.append(df_s)
        else:
            log.warning("ssGSEA score file missing for %s — skipped", cohort)

    if ssgsea_parts:
        combined_ssgsea = pd.concat(ssgsea_parts)
        combined_ssgsea.index.name = "sample_id"
        combined_ssgsea_path = results_dir / "combined_ssgsea_scores.tsv"
        combined_ssgsea.to_csv(combined_ssgsea_path, sep="\t")
        log.info("Combined ssGSEA scores saved → %s", combined_ssgsea_path)
    else:
        log.warning("No ssGSEA score files found — combined_ssgsea_scores.tsv not created")

    # ── combine and save ──────────────────────────────────────────────────────
    if all_labels:
        combined = pd.concat(all_labels, ignore_index=True)
        combined_path = results_dir / "combined_classification.tsv"
        combined.to_csv(combined_path, sep="\t", index=False)
        log.info("=" * 60)
        log.info("Combined classification saved → %s", combined_path)
        log.info("Total samples classified: %d across %d cohort(s)",
                 len(combined), combined["cohort"].nunique())
        log.info("TME distribution (all cohorts):\n%s",
                 combined["TME"].value_counts().to_string())

    if failed:
        log.warning("Failed cohorts (%d): %s", len(failed), failed)
    else:
        log.info("All cohorts completed successfully.")


if __name__ == "__main__":
    main()
