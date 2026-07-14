# Methylation Trio Platform

Offline Streamlit application for exploratory nanopore methylation DMR analysis
of one proband and two relatives. Roles, sexes, region scope, and every analysis
threshold are configurable. Input files are selected from a mounted data
directory rather than uploaded, and reference assemblies are managed for you.

## Start with Docker

Requirements: Docker Compose and a data directory containing the three modBAMs
(with `.bai` indexes).

```bash
export METHYL_TRIO_DATA=/absolute/path/to/data
docker compose up --build
```

Open <http://localhost:8501>. The host data directory is mounted read-only at
`/data`; results are written under `./runs`; downloaded references persist in the
`methyl-trio-references` volume mounted at `/references`.

The image includes modkit, samtools/htslib, methylartist, modbamtools, and
Python. The reference FASTA, GENCODE annotation, and UCSC CpG islands for the
selected assembly (hg38 or hg19) are downloaded and cached once on first use,
then reused offline. Set `METHYL_TRIO_REFERENCE_CACHE` to relocate the cache.

## Workflow and method

1. Setup: choose the three BAMs from the data directory with a rooted file
   browser (no uploads), pick the genome assembly (hg38/hg19) and prepare it
   once, optionally select a phased VCF, and validate BAM/reference contig
   lengths, basecaller models, and HP-tag availability. Optional relationship and
   clinical-status fields include inline explanations.
2. Regions supports a genome-wide CpG-island scan, selected chromosomes, or
   GENCODE-derived promoter/gene-body intervals for an editable gene panel.
3. Thresholds renders controls and rationale from `core/thresholds.py`.
4. Run generates indexed bedMethyl with `modkit pileup`, scores the same regions
   for P-vs-R1, P-vs-R2, and R1-vs-R2 with `modkit dmr pair`, and ranks results.
5. Results provides a table, plot, verdict, caveats, a self-contained HTML
   report, and a complete-run ZIP archive.

The empirical null is the configured percentile (99th by default) of absolute
R1-vs-R2 effect sizes. A candidate must:

- overlap in both proband comparisons;
- have the same effect direction;
- exceed the null cutoff in both comparisons;
- pass site-count and p-value thresholds; and
- not overlap a region observed in the relative-relative null.

When exactly one relative is affected and one is unaffected, the app switches
to phenotype-segregation ranking: affected samples must be similar to each
other and both must differ concordantly from the unaffected relative. Two
explicitly affected relatives cannot support this analysis because there is no
unaffected comparator. Missing clinical status remains optional and falls back
to the original exploratory trio design with a report caveat.

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
# The picker is rooted at METHYL_TRIO_DATA_ROOT (defaults to your home directory
# when unset); references cache under METHYL_TRIO_REFERENCE_CACHE.
METHYL_TRIO_DATA_ROOT=/path/to/data \
  PYTHONPATH="$PWD" streamlit run app/streamlit_app.py
```

## Outputs

Each run writes:

- `config.json` — resolved run manifest;
- `pipeline.log` — tool output;
- three pairwise DMR tables and sample bedMethyl files;
- `proband_specific_DMRs.tsv`;
- `dmr_effects.png`, `summary.json`, and `report.html`.

The Results page can also build `methyl_trio_run.zip`, a deterministic archive of
the entire run directory (manifest, logs, pileups, pairwise tables, ranked
candidates, figures, and report) for download.

## Scientific limitations

This is a family prioritization screen, not a cohort analysis or diagnostic
test. Three related samples cannot estimate population variance. Basecalling,
sequencing, and processing batches can mimic methylation effects. Missing
parental lineages leave parent-of-origin mQTL effects unresolved; the report
generates a composition-specific warning. HP tags are optional for DMR calling
but required for haplotype read-level interpretation.

A phased VCF and identified mother and father make parent-of-origin and mQTL
follow-up possible, but their presence alone does not establish either effect.
The report labels each interpretive question (phenotype segregation, parent of
origin, mQTL) by whether the required metadata is available, rather than treating
absent inputs as evidence of no effect.

## Validation status

Unit and synthetic smoke tests cover trio validation, sex-chromosome rules,
null-threshold intersection/ranking, `N_other` correction, command construction,
threshold bounds, targeted statistics, and HTML export.

BH11998 parity is not claimed by this repository: the original `00`–`08`
scripts, external-drive BAMs, and expected `proband_specific_DMRs.tsv` were not
present in the source repository. To perform the requested parity check, mount
those inputs, run the pilot chromosomes (`chr1`, `chr2`, `chr11`, `chr15`), and
compare coordinates/effects against the historical output.
