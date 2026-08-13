"""Build the frozen KOLF-common LINCS payload from the original Level-5 GCTX.

This is a release-data build script, not runtime package code.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


def decode(values):
    return np.asarray([
        value.decode("utf-8") if isinstance(value, bytes) else str(value)
        for value in values
    ])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gctx", required=True)
    parser.add_argument("--sig-info", required=True)
    parser.add_argument("--basis", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--chunk-size", type=int, default=3000)
    args = parser.parse_args()

    basis = pd.read_csv(args.basis, sep="\t", index_col=0)
    basis.index = basis.index.astype(str)
    weights = basis.to_numpy(dtype=np.float64)
    norms = np.sqrt(np.sum(weights * weights, axis=0))
    if np.any(~np.isfinite(norms)) or np.any(norms <= 0):
        raise ValueError("Invalid KOLF basis norm")
    weights /= norms

    metadata = pd.read_csv(args.sig_info, sep="\t", low_memory=False)
    required = {
        "sig_id", "pert_id", "pert_iname", "pert_type", "cell_id",
        "pert_idose", "pert_itime",
    }
    missing = required.difference(metadata.columns)
    if missing:
        raise ValueError(f"Missing LINCS metadata columns: {sorted(missing)}")
    metadata = metadata.loc[
        (metadata["pert_type"].astype(str) == "trt_cp")
        & (metadata["pert_itime"].astype(str) == "24 h")
    ].copy()

    with h5py.File(args.gctx, "r") as handle:
        matrix = handle["0/DATA/0/matrix"]
        sig_ids = decode(handle["0/META/COL/id"][:])
        gene_ids = decode(handle["0/META/ROW/id"][:])
        sig_to_index = {value: i for i, value in enumerate(sig_ids)}
        gene_to_index = {value: i for i, value in enumerate(gene_ids)}

        missing_genes = basis.index.difference(gene_to_index)
        if len(missing_genes):
            raise ValueError(f"{len(missing_genes)} basis genes absent from GCTX")
        gene_pairs = sorted(
            ((gene_to_index[gene], row) for row, gene in enumerate(basis.index)),
            key=lambda item: item[0],
        )
        gctx_gene_indices = np.asarray([item[0] for item in gene_pairs])
        basis_row_indices = np.asarray([item[1] for item in gene_pairs])
        ordered_weights = weights[basis_row_indices, :]

        metadata = metadata.loc[metadata["sig_id"].isin(sig_to_index)].copy()
        metadata["gctx_index"] = metadata["sig_id"].map(sig_to_index).astype(int)
        metadata.sort_values("gctx_index", inplace=True)
        projected_blocks = []
        total = len(metadata)
        for start in range(0, total, args.chunk_size):
            stop = min(start + args.chunk_size, total)
            block_meta = metadata.iloc[start:stop]
            indices = block_meta["gctx_index"].to_numpy()
            expression = matrix[indices, :][:, gctx_gene_indices]
            projected = np.asarray(expression, dtype=np.float64) @ ordered_weights
            block = pd.DataFrame(projected, columns=basis.columns)
            block.insert(0, "sig_id", block_meta["sig_id"].to_numpy())
            projected_blocks.append(block)
            print(f"Projected {stop}/{total}", flush=True)

    projected = pd.concat(projected_blocks, ignore_index=True)
    instance = metadata.drop(columns=["gctx_index"]).merge(
        projected, on="sig_id", how="inner", validate="one_to_one"
    )
    group = ["pert_id", "pert_iname", "cell_id"]
    scores = instance.groupby(group, sort=False)[list(basis.columns)].median()
    counts = instance.groupby(group, sort=False).agg(
        n_signatures=("sig_id", "size"),
        n_doses=("pert_idose", "nunique"),
        n_times=("pert_itime", "nunique"),
        pert_times=("pert_itime", lambda x: "|".join(sorted(set(map(str, x))))),
    )
    result = counts.join(scores).reset_index().rename(
        columns={"cell_id": "cell_iname"}
    )
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, sep="\t", index=False, compression="gzip")
    print(
        f"Wrote {len(result)} drug-cell profiles and {len(basis.columns)} programs",
        flush=True,
    )


if __name__ == "__main__":
    main()
