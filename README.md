# Mechanism-Aware Penetrance Prediction

A tool that predicts a gene's (and its variants') **tendency toward incomplete
vs. full penetrance** from **mechanism + constraint + the pathogenic
allele-frequency spectrum**, and uses that learned gene-level prediction as the
**Beta prior** for an off-the-shelf per-variant Bayesian layer.

The center of gravity is the part the literature flags as unsolved: constraint
metrics (pLI, LOEUF, `shet`) measure selection against heterozygous loss of
function, which is repeatedly warned to be **not** penetrance (Fuller et al.).
No published model predicts, from a gene, *how incompletely penetrant its
pathogenic variants tend to be*, grounded in mechanism. That is the
tubulin-vs-BRCA question:

- **Tubulins** (`TUBA1A`, `TUBB`, `TUBB3`): pathogenic variants act by
  **dominant-negative / gain-of-function** on the microtubule lattice → large,
  near-deterministic effect → **high penetrance**.
- **Haploinsufficiency / dosage genes** (`BRCA1/2`, `LDLR`, `SDHx`, sarcomere,
  channelopathies): a single hit leaves a functioning allele → effect modifiable
  by background, environment and age → **incomplete / age-dependent penetrance**.

## What is novel vs. reused

| Layer | Status | Source |
|---|---|---|
| Per-variant Beta-Binomial posterior | **reused** | CalPen / Kroncke, PLoS Genet 2020 |
| Case/control + prevalence Bayes | **reused** | ADpenetrance (Wright/KCL) |
| Max credible allele-frequency bound | **reused** | Whiffin/Ware, Genet Med 2017 |
| LLM literature mining of carrier counts | **reused** (optional adapter) | GeneVariantFetcher |
| **Learned, mechanism-aware gene penetrance-propensity model that supplies the Beta prior** | **new** | this repo |

## Architecture

```
Gene-level inputs                         Per-variant layer (off the shelf)
  mechanism (LoF/HI vs DN/GoF, inferred)     carrier counts (ClinVar/gnomAD
  constraint (LOEUF, pLI, shet, pHaplo,        or optional literature adapter)
    pTriplo, missense o/e)                    Whiffin/Ware max-credible-AF bound
  pathogenic allele-frequency spectrum
             │                                        │
             ▼                                        ▼
   Gene incomplete-penetrance   ──►  Beta prior  ──►  Beta-Binomial posterior
   propensity model (GBM)            (alpha, beta)          │
                                                            ▼
                             penetrance: point estimate + credible interval + tier
```

## Package layout

```
src/penetrance/
  labels/       Component 1 - curated literature-derived penetrance ground truth
  features/     Component 2 - mechanism inference + mechanism-aware feature matrix
  gene_model/   Component 3 - GBM regression, family-aware CV, calibration, Beta-prior mapping  (core contribution)
  variant/      Component 4 - Beta-Binomial posterior, case/control Bayes, Whiffin/Ware bound
  adapters/     Component 5 - pluggable carrier-count sources (ClinVar/gnomAD + literature)
  pipeline.py   end-to-end orchestration (gene prior -> per-variant posterior)
  validation.py Component 6 - the experiments that prove the point
  cli.py        command-line interface
data/ (packaged) genes.csv, variants.csv
tests/          unit + integration tests
```

## Install

```bash
pip install -e .           # core (numpy/pandas/scipy/scikit-learn)
pip install -e ".[dev]"    # adds lightgbm + shap + pytest
```

LightGBM is used for the gene model when available; otherwise the package falls
back to scikit-learn's `HistGradientBoostingRegressor` automatically.

## Usage

```bash
# Train the gene model and report gene-family-aware cross-validation metrics
penetrance train

# Predict a gene's penetrance propensity and the Beta prior it induces
penetrance gene TUBA1A BRCA1 LDLR HFE

# Estimate a variant's penetrance with the mechanism prior vs a flat prior
penetrance variant HFE:p.Cys282Tyr SDHB:p.Ile127Ser --provenance

# Run the validation experiments
penetrance validate
```

Programmatically:

```python
from penetrance.pipeline import PenetrancePipeline

pipe = PenetrancePipeline().fit()
pipe.predict_gene("BRCA1")                       # propensity + Beta prior
est = pipe.estimate_variant("SDHB:p.Ile127Ser")  # posterior + credible interval + tier
print(est.point_estimate, (est.ci_low, est.ci_high), est.confidence_tier)
```

## Validation results (packaged curated label set)

Running `penetrance validate` on the shipped labels:

- **Gene model (gene-family-aware CV, out-of-fold):** Spearman ≈ 0.74,
  Pearson ≈ 0.67, MAE ≈ 0.16, expected calibration error ≈ 0.05. The dominant
  drivers (both split-gain importance and mean-|SHAP|) are mechanism (inferred
  `dn_gof_score` and the missense:LoF pathogenic ratio), constraint (pLI,
  pTriplo, LOEUF) and the pathogenic allele-frequency spectrum — i.e. the model
  learns *mechanism + constraint + AF spectrum*, not gene identity (paralogs are
  held out together).
- **Mechanism separation:** mean predicted propensity ≈ 0.93 for high-penetrance
  DN/GoF families (tubulin/collagen/fibrillin) vs ≈ 0.41 for incomplete-penetrance
  families (HBOC/FH/paraganglioma) — a ≈ 0.52 gap.
- **Sparse-variant prior benefit** (mask-and-recover): the mechanism prior cuts
  the mean absolute error of the penetrance estimate vs. a flat `Beta(1,1)` prior
  by ≈ **54%** with a single observed carrier, tapering to ≈ 28% at 8 carriers —
  exactly where a mechanism-informed prior should help and then get out of the way.

## Data provenance and honesty

The packaged `data/genes.csv` and `data/variants.csv` are **curated
approximations** drawn from the population/biobank and clinical literature cited
in each row's `source` field (Wright et al. Nat Genet 2024; Forrest/Huang et al.
Nat Genet 2025; the Science 2024/2025 ML-penetrance gene set; ClinGen/OMIM; and
the per-gene studies named). They make the pipeline runnable and reproducible
offline. The adapters (`penetrance.adapters`) are the seam for wiring in live
gnomAD/ClinVar frequencies and literature-mined counts:

- `FrequencyCountAdapter` — point it at a table pulled from live ClinVar + gnomAD.
- `LiteratureCountAdapter` — a fork-of-GeneVariantFetcher miner with a **pluggable
  extractor** (drop in an LLM), a **variant-matching gate**, **ontology
  normalization**, and **provenance** quotes. Ships with a dependency-free
  regex extractor so it runs offline.

## Key references

- Wright et al., *Guidance for estimating penetrance ... population cohorts*, Nat Genet 2024.
- Forrest, Huang et al., *Using large-scale population-based data ...*, Nat Genet 2025.
- *Machine learning-based penetrance of genetic variants*, Science 2024/2025.
- CalPen (PLOS); ADpenetrance (Wright/KCL); Kroncke Bayesian, PLoS Genet 2020.
- Fuller et al., *Measuring intolerance to mutation*; GeneBayes `shet` (constraint ≠ penetrance).
- Whiffin/Ware, *maximum credible allele frequency*, Genet Med 2017.
- GeneVariantFetcher (kroncke-lab) — optional literature-count adapter only.
