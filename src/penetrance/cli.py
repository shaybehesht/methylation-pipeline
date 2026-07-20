"""Command-line interface for the penetrance tool."""

from __future__ import annotations

import argparse
import sys
from typing import List, Optional

import pandas as pd

from penetrance.pipeline import PenetrancePipeline
from penetrance.validation import (
    gene_model_calibration,
    mechanism_separation,
    sparse_prior_benefit,
)


def _build_pipeline() -> PenetrancePipeline:
    pipe = PenetrancePipeline()
    pipe.fit()
    return pipe


def _cmd_train(args) -> int:
    pipe = _build_pipeline()
    cv = pipe.cv_result_
    print(f"Gene model backend: {pipe.model.backend_}")
    print(f"Trained on {len(pipe.labels.genes)} genes across "
          f"{pipe.labels.genes['gene_family'].nunique()} families.")
    print("\nFamily-aware CV metrics (out-of-fold):")
    for k, v in cv.metrics.items():
        print(f"  {k:10s}: {v:.4f}")
    print("\nTop feature importances:")
    print(pipe.model.feature_importance().head(8).to_string())
    return 0


def _cmd_gene(args) -> int:
    pipe = _build_pipeline()
    for gene in args.gene:
        pred = pipe.predict_gene(gene)
        if pred is None:
            print(f"{gene}: not found in label set")
            continue
        prior = pred.prior
        print(
            f"{gene}: propensity={pred.propensity:.3f} (+/-{pred.std:.3f})  "
            f"Beta prior=({prior.alpha:.2f}, {prior.beta:.2f}), "
            f"prior mean={prior.mean:.3f}, strength={prior.strength:.1f}"
        )
    return 0


def _cmd_variant(args) -> int:
    pipe = _build_pipeline()
    for vid in args.variant_id:
        counts = pipe.carrier_counts(vid)
        mech = pipe.estimate_variant(vid, use_gene_prior=True)
        flat = pipe.estimate_variant(vid, use_gene_prior=False)
        if mech is None:
            print(f"{vid}: not found")
            continue
        print(f"\n=== {vid} ===")
        if counts is not None:
            print(f"  carriers: affected={counts.affected:g}, unaffected={counts.unaffected:g}, "
                  f"source={counts.source}")
        print(f"  mechanism prior: p={mech.point_estimate:.3f} "
              f"[{mech.ci_low:.3f}, {mech.ci_high:.3f}] tier={mech.confidence_tier}")
        print(f"  flat prior     : p={flat.point_estimate:.3f} "
              f"[{flat.ci_low:.3f}, {flat.ci_high:.3f}] tier={flat.confidence_tier}")
        if mech.af_upper_bound is not None:
            print(f"  Whiffin/Ware max-credible-AF penetrance bound: {mech.af_upper_bound:.3f}")
        if args.provenance and counts is not None:
            for p in counts.provenance:
                print(f"    - {p}")
    return 0


def _cmd_validate(args) -> int:
    pipe = _build_pipeline()
    pd.set_option("display.width", 120)

    print("== Gene-model CV calibration ==")
    calib = gene_model_calibration(pipe)
    for k, v in calib["metrics"].items():
        print(f"  {k:10s}: {v:.4f}")
    print(f"  expected calibration error: {calib['expected_calibration_error']:.4f}")
    print(calib["reliability"].to_string(index=False))

    print("\n== Mechanism separation (tubulin/collagen/fibrillin vs hboc/fh/pgl) ==")
    sep = mechanism_separation(pipe)
    for k, v in sep.items():
        print(f"  {k}: {v:.3f}")

    print("\n== Sparse-variant prior benefit (flat vs mechanism) ==")
    res = sparse_prior_benefit(pipe, n_trials=args.trials)
    summary = res.summary.copy()
    summary["mae_reduction"] = res.improvement().values
    print(summary.to_string(index=False))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="penetrance",
        description="Mechanism-aware penetrance prediction (gene prior + per-variant Bayesian layer).",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p_train = sub.add_parser("train", help="Train the gene model and report CV metrics.")
    p_train.set_defaults(func=_cmd_train)

    p_gene = sub.add_parser("gene", help="Predict a gene's penetrance propensity + Beta prior.")
    p_gene.add_argument("gene", nargs="+", help="Gene symbol(s).")
    p_gene.set_defaults(func=_cmd_gene)

    p_var = sub.add_parser("variant", help="Estimate penetrance for a variant.")
    p_var.add_argument("variant_id", nargs="+", help="Variant id(s) from the label set.")
    p_var.add_argument("--provenance", action="store_true", help="Print carrier-count provenance.")
    p_var.set_defaults(func=_cmd_variant)

    p_val = sub.add_parser("validate", help="Run the validation experiments.")
    p_val.add_argument("--trials", type=int, default=300, help="Monte-Carlo trials for the sparse experiment.")
    p_val.set_defaults(func=_cmd_validate)

    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv if argv is not None else sys.argv[1:])
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
