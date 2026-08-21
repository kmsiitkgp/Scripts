#!/usr/bin/env python3
"""
SCRNA_run_velocity.py

Runs RNA velocity on one lineage's velocity_input.h5ad using BOTH scVelo
(dynamical model) and veloVI (scvi-tools), sharing one preprocessing pass
(filter/normalize + moments) since both methods need it.

Design: each method's block is wrapped independently and writes its own
output file only if it succeeds. Either method can be commented out below
(or can fail/raise) without affecting the other — the .nf process's outputs
are declared `optional: true` for exactly this reason, so no Nextflow-side
change is needed if you comment out or disable one method.

  scVelo  -> {label}_scvelo.h5ad + {label}_scvelo_stream.png
  veloVI  -> {label}_velovi.h5ad + {label}_velovi_stream.png

Usage:
    python SCRNA_run_velocity.py --input velocity_input.h5ad --label Epithelial --output-dir .
"""

import argparse
import sys
import traceback

import matplotlib
matplotlib.use("Agg")  # headless — no display available in a container/HPC job
import matplotlib.pyplot as plt
import scvelo as scv


def parse_args():
    parser = argparse.ArgumentParser(description="Run scVelo + veloVI RNA velocity on one lineage")
    parser.add_argument("--input", required=True, help="Path to velocity_input.h5ad")
    parser.add_argument("--label", required=True, help="Lineage label, e.g. Epithelial (used in plot titles/filenames)")
    parser.add_argument("--output-dir", required=True, help="Directory to write outputs into")
    parser.add_argument("--n-top-genes", type=int, default=2000, help="Highly variable genes to keep (default: 2000)")
    parser.add_argument("--n-pcs", type=int, default=30, help="PCs for neighbor graph (default: 30)")
    parser.add_argument("--n-neighbors", type=int, default=30, help="Neighbors for moments (default: 30)")
    parser.add_argument("--velovi-max-epochs", type=int, default=500, help="Max training epochs for veloVI (default: 500)")
    return parser.parse_args()


def load_and_preprocess(input_path, n_top_genes, n_pcs, n_neighbors):
    """Shared preprocessing used by both scVelo and veloVI."""
    print(f"Reading: {input_path}")
    adata = scv.read(input_path)

    required_layers = {"mature", "nascent", "ambiguous"}
    missing = required_layers - set(adata.layers.keys())
    if missing:
        sys.exit(
            f"❌ ERROR: expected layers {required_layers} from kb-python nac workflow, "
            f"missing: {missing}. Found layers: {list(adata.layers.keys())}"
        )

    # scVelo/veloVI both expect "spliced"/"unspliced" layer names — kb-python's
    # nac workflow names them "mature"/"nascent"/"ambiguous". Map explicitly;
    # fold ambiguous into spliced per kb-python's documented nac -> scVelo
    # handoff recommendation (ambiguous reads are usually mostly exonic).
    adata.layers["spliced"] = adata.layers["mature"]
    adata.layers["unspliced"] = adata.layers["nascent"]
    if "ambiguous" in adata.layers:
        adata.layers["spliced"] = adata.layers["spliced"] + adata.layers["ambiguous"]

    print(f"Input: {adata.n_obs} cells x {adata.n_vars} genes")
    frac_s = adata.layers["spliced"].sum() / (adata.layers["spliced"].sum() + adata.layers["unspliced"].sum())
    print(f"Spliced fraction: {frac_s:.2%} (sanity check — typically ~0.7-0.95 for nac-derived counts)")

    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=n_top_genes)
    scv.pp.moments(adata, n_pcs=n_pcs, n_neighbors=n_neighbors)

    return adata


def run_scvelo(adata, label, output_dir):
    """scVelo dynamical model. Writes {label}_scvelo.h5ad + stream plot on success."""
    print("\n=== Running scVelo (dynamical model) ===")
    adata = adata.copy()  # isolate from veloVI's mutations of the shared object

    scv.tl.recover_dynamics(adata, n_jobs=-1)
    scv.tl.velocity(adata, mode="dynamical")
    scv.tl.velocity_graph(adata)
    scv.tl.latent_time(adata)

    if "X_umap" in adata.obsm:
        scv.pl.velocity_embedding_stream(
            adata, basis="umap", color="CellType",
            title=f"{label} — scVelo (dynamical)",
            save=f"{label}_scvelo_stream.png", show=False,
        )
        plt.close("all")
    else:
        print("⚠️  WARNING: no X_umap in obsm — skipping scVelo embedding plot.", file=sys.stderr)

    out_path = f"{output_dir}/{label}_scvelo.h5ad"
    adata.write_h5ad(out_path)
    print(f"✅ scVelo done: wrote {out_path} ({adata.n_obs} cells x {adata.n_vars} genes)")


def run_velovi(adata, label, output_dir, max_epochs):
    """veloVI deep generative velocity model. Writes {label}_velovi.h5ad + stream plot on success."""
    print("\n=== Running veloVI ===")
    import scvi
    from velovi import VELOVI

    adata = adata.copy()  # isolate from scVelo's mutations of the shared object

    VELOVI.setup_anndata(adata, spliced_layer="Ms", unspliced_layer="Mu")
    model = VELOVI(adata)
    model.train(max_epochs=max_epochs)

    adata.layers["velocity"] = model.get_velocity()
    adata.obs["latent_time_velovi"] = model.get_latent_time().mean(axis=1)

    scv.tl.velocity_graph(adata, vkey="velocity")

    if "X_umap" in adata.obsm:
        scv.pl.velocity_embedding_stream(
            adata, basis="umap", color="CellType", vkey="velocity",
            title=f"{label} — veloVI",
            save=f"{label}_velovi_stream.png", show=False,
        )
        plt.close("all")
    else:
        print("⚠️  WARNING: no X_umap in obsm — skipping veloVI embedding plot.", file=sys.stderr)

    out_path = f"{output_dir}/{label}_velovi.h5ad"
    adata.write_h5ad(out_path)
    print(f"✅ veloVI done: wrote {out_path} ({adata.n_obs} cells x {adata.n_vars} genes)")


def main():
    args = parse_args()
    scv.settings.verbosity = 3
    scv.settings.figdir = args.output_dir

    adata = load_and_preprocess(args.input, args.n_top_genes, args.n_pcs, args.n_neighbors)

    # -------------------------------------------------------------------
    # Each method is independently try/excepted: a veloVI failure (or it
    # being commented out entirely) never blocks scVelo's already-written
    # output, and vice versa. The .nf process declares both outputs
    # `optional: true` for exactly this reason — no Nextflow change needed
    # if you comment out either block below.
    # -------------------------------------------------------------------
    try:
        run_scvelo(adata, args.label, args.output_dir)
    except Exception:
        print(f"❌ ERROR: scVelo failed for {args.label}", file=sys.stderr)
        traceback.print_exc()

    try:
        run_velovi(adata, args.label, args.output_dir, args.velovi_max_epochs)
    except Exception:
        print(f"❌ ERROR: veloVI failed for {args.label}", file=sys.stderr)
        traceback.print_exc()


if __name__ == "__main__":
    main()
