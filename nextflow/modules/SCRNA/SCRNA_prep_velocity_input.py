#!/usr/bin/env python3
"""
SCRNA_prep_velocity_input.py

Subsets the merged kb-python AnnData to one lineage's barcodes (as determined
by that lineage's own Seurat _Final object), and attaches CellType/cluster/UMAP
columns from barcodes.csv (written by SCRNA_annotate_clusters.R alongside
annotated_seurat.rds).

kb-python's --filter bustools cell-calling and Cell Ranger's cell-calling are
independent algorithms — barcode sets will not match perfectly in either
direction. This performs a strict inner join and reports the overlap so a
badly mismatched sample/chemistry setup is caught early rather than silently
producing an undersized velocity input.

Usage:
    python SCRNA_prep_velocity_input.py --merged merged.h5ad --barcodes barcodes.csv --output velocity_input.h5ad
"""

import argparse
import sys

import anndata as ad
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Subset merged h5ad to one lineage's barcodes")
    parser.add_argument("--merged", required=True, help="Path to merged.h5ad (all samples)")
    parser.add_argument("--barcodes", required=True, help="Path to barcodes.csv (from SCRNA_annotate_clusters.R)")
    parser.add_argument("--output", required=True, help="Path to write lineage velocity_input.h5ad")
    return parser.parse_args()


def main():
    args = parse_args()

    print(f"Reading merged AnnData: {args.merged}")
    merged = ad.read_h5ad(args.merged)

    print(f"Reading lineage barcodes: {args.barcodes}")
    barcodes_df = pd.read_csv(args.barcodes, dtype={"barcode": str})
    barcodes_df = barcodes_df.set_index("barcode")

    seurat_barcodes = set(barcodes_df.index)
    merged_barcodes = set(merged.obs_names)

    common = seurat_barcodes & merged_barcodes
    only_seurat = seurat_barcodes - merged_barcodes

    overlap_pct = 100 * len(common) / len(seurat_barcodes) if seurat_barcodes else 0

    print(
        f"Barcode overlap: {len(common)} / {len(seurat_barcodes)} Seurat barcodes "
        f"found in merged h5ad ({overlap_pct:.1f}%)"
    )
    print(f"  Seurat barcodes missing from h5ad (kb-python filter excluded them): {len(only_seurat)}")

    if len(common) == 0:
        sys.exit(
            "❌ ERROR: zero barcode overlap between Seurat object and merged h5ad — "
            "check sample_id naming consistency between Cell Ranger and kb-python runs"
        )

    if overlap_pct < 80:
        print(
            f"⚠️  WARNING: barcode overlap ({overlap_pct:.1f}%) is lower than expected (<80%). "
            "This may indicate a sample_id mismatch or chemistry misconfiguration — "
            "verify before trusting downstream velocity results.",
            file=sys.stderr,
        )

    # Strict inner join: subset merged AnnData to only the common, ordered
    # barcodes, then attach this lineage's CellType/cluster/UMAP columns in
    # that same order.
    common_ordered = [bc for bc in merged.obs_names if bc in common]
    subset = merged[common_ordered].copy()

    barcodes_aligned = barcodes_df.loc[common_ordered]
    subset.obs["CellType"] = barcodes_aligned["CellType"].values
    subset.obs["cluster"] = barcodes_aligned["cluster"].astype(str).values
    if {"UMAP_1", "UMAP_2"}.issubset(barcodes_aligned.columns):
        subset.obsm["X_umap"] = barcodes_aligned[["UMAP_1", "UMAP_2"]].values
    else:
        print("⚠️  WARNING: barcodes.csv has no UMAP_1/UMAP_2 columns — X_umap not set.", file=sys.stderr)

    print(f"Lineage velocity input: {subset.n_obs} cells x {subset.n_vars} genes")
    subset.write_h5ad(args.output)
    print(f"✅ SUCCESS: wrote {args.output}")


if __name__ == "__main__":
    main()
