# Methylation Trio Platform

Offline Streamlit application for exploratory nanopore methylation DMR analysis
of one proband and two relatives. Roles, sexes, paths, region scope, and every
analysis threshold are configurable.

## Start with Docker

Requirements: Docker Compose and a reference FASTA matching all three BAM
headers.

```bash
export METHYL_TRIO_DATA=/absolute/path/to/data
docker compose up --build
```

Open <http://localhost:8501>. The host data directory is mounted read-only at
`/data`; results are written under `./runs`.

The image includes modkit, samtools/htslib, methylartist, modbamtools, Python,
GENCODE annotation, and UCSC hg38 CpG islands. Annotation downloads occur only
during image build. Runtime analysis does not require internet access.

## Workflow and method

1. Setup validates exactly three unique samples, one proband, BAM/reference
   contig lengths, reported basecaller models, and HP-tag availability.
2. Regions supports a genome-wide CpG-island scan, selected chromosomes, or
   GENCODE-derived promoter/gene-body intervals for an editable gene panel.
3. Thresholds renders controls and rationale from `core/thresholds.py`.
4. Run generates indexed bedMethyl with `modkit pileup`, scores the same regions
   for P-vs-R1, P-vs-R2, and R1-vs-R2 with `modkit dmr pair`, and ranks results.
5. Results provides a table, plot, verdict, caveats, and self-contained HTML.

The empirical null is the configured percentile (99th by default) of absolute
R1-vs-R2 effect sizes. A candidate must:

- overlap in both proband comparisons;
- have the same effect direction;
- exceed the null cutoff in both comparisons;
- pass site-count and p-value thresholds; and
- not overlap a region observed in the relative-relative null.

Region p-values are two-sided Fisher exact tests reconstructed from modkit
modified and unmodified totals. Targeted runs additionally apply the configured
effect-size and alpha cutoffs. `core/analysis.targeted_test` exposes a
Mann–Whitney/Wilcoxon-rank-sum test for per-site targeted summaries.

chrX is valid only when both samples in a comparison are female, chrY only when
both are male, and mitochondrial contigs are excluded. These rules are applied
independently to all three comparisons.

## Local development

```bash
python -m pip install -e '.[test]'
pytest
streamlit run app/streamlit_app.py
```

## Outputs

Each run writes:

- `config.json` — resolved run manifest;
- `pipeline.log` — tool output;
- three pairwise DMR tables and sample bedMethyl files;
- `proband_specific_DMRs.tsv`;
- `dmr_effects.png`, `summary.json`, and `report.html`.

## Scientific limitations

This is a family prioritization screen, not a cohort analysis or diagnostic
test. Three related samples cannot estimate population variance. Basecalling,
sequencing, and processing batches can mimic methylation effects. Missing
parental lineages leave parent-of-origin mQTL effects unresolved; the report
generates a composition-specific warning. HP tags are optional for DMR calling
but required for haplotype read-level interpretation.

## Validation status

Unit and synthetic smoke tests cover trio validation, sex-chromosome rules,
null-threshold intersection/ranking, `N_other` correction, command construction,
threshold bounds, targeted statistics, and HTML export.

BH11998 parity is not claimed by this repository: the original `00`–`08`
scripts, external-drive BAMs, and expected `proband_specific_DMRs.tsv` were not
present in the source repository. To perform the requested parity check, mount
those inputs, run the pilot chromosomes (`chr1`, `chr2`, `chr11`, `chr15`), and
compare coordinates/effects against the historical output.
