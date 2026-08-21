#!/usr/bin/env python3
"""
merge_h5ad.py

Merges per-sample kb-python (--h5ad, --workflow nac) output into a single
AnnData, and rewrites barcodes to match the Seurat object's convention:

    kb-python barcode:  AAACCAAAGAAGCCTC
    Seurat barcode:      {sample_id}_AAACCAAAGAAGCCTC-1

(Seurat's format comes from: Read10X(strip.suffix=FALSE) keeps Cell Ranger's
"-1" suffix -> merge(add.cell.ids = sample_ids) prepends "{sample_id}_".
See SCRNA_load_cellranger_and_classify_droplets.R / SCRNA_merge_and_plot_qc.R.)

Usage:
    python merge_h5ad.py --output merged.h5ad SAMPLE1=sample1.h5ad SAMPLE2=sample2.h5ad ...
"""

import argparse
import sys

import anndata as ad


def parse_args():
    parser = argparse.ArgumentParser(description="Merge per-sample kb-python h5ad files")
    parser.add_argument("--output", required=True, help="Path to write merged .h5ad")
    parser.add_argument(
        "sample_files",
        nargs="+",
        help="SAMPLE_ID=path.h5ad pairs, one per sample",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    adatas = {}
    for entry in args.sample_files:
        if "=" not in entry:
            sys.exit(f"❌ ERROR: expected SAMPLE_ID=path.h5ad, got '{entry}'")
        sample_id, path = entry.split("=", 1)

        print(f"Reading {sample_id}: {path}")
        adata = ad.read_h5ad(path)

        # Harmonize barcodes to match the Seurat object's convention:
        # kb-python barcode "AAACCAAAGAAGCCTC" -> "{sample_id}_AAACCAAAGAAGCCTC-1"
        adata.obs_names = [f"{sample_id}_{bc}-1" for bc in adata.obs_names]
        adata.obs["Sample"] = sample_id

        adatas[sample_id] = adata

    print(f"Concatenating {len(adatas)} samples...")
    merged = ad.concat(
        adatas,
        join="outer",       # union of genes across samples; missing genes filled 0
        index_unique=None,  # barcodes are already uniquely prefixed per-sample above
    )

    # obs_names should now be unique across the whole merged object; verify rather
    # than silently letting a downstream barcode-subset step fail confusingly later
    n_total = merged.n_obs
    n_unique = merged.obs_names.nunique()
    if n_total != n_unique:
        sys.exit(
            f"❌ ERROR: {n_total - n_unique} duplicate barcodes after merge "
            f"({n_total} total, {n_unique} unique) — check sample_id uniqueness"
        )

    print(f"Merged AnnData: {merged.n_obs} cells x {merged.n_vars} genes")
    merged.write_h5ad(args.output)
    print(f"✅ SUCCESS: wrote {args.output}")


if __name__ == "__main__":
    main()
