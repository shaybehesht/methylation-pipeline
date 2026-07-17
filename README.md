# 🥭 MANGO — Methylation Analysis for Novel Genomic Outcomes

MANGO (Methylation Analysis for Novel Genomic Outcomes) is an offline Streamlit
application for exploratory nanopore methylation DMR analysis of one proband and
two relatives. Roles, sexes, region scope, and every analysis threshold are
configurable. Local modBAMs are selected through a server-side file browser, and
reference assemblies are downloaded once and cached for offline reuse.

## Start with Docker

Requirements: Docker Compose. No reference FASTA needs to be supplied by hand —
choose hg38 or hg19 in the app and it is fetched and prepared into a persistent
cache the first time it is used.

```bash
export METHYL_TRIO_DATA=/absolute/path/to/data
export METHYL_TRIO_REFERENCE_CACHE=/absolute/path/to/reference-cache
docker compose up --build
```

Open <http://localhost:8501>. The host data directory is mounted read-only at
`/data` (browse it in Setup); the reference cache is mounted read-write at
`/reference-cache`; results are written under `./runs`.

The image includes modkit, samtools/htslib, methylartist, modbamtools, and
Python. Reference FASTAs and matching GENCODE/CpG-island annotations are **not**
baked into the image; they are downloaded on demand into the mounted reference
cache so the image stays small and both assemblies are not duplicated. Only the
one-time reference preparation needs internet access — subsequent analysis runs
are fully offline as long as the cache is persisted.

## Install and share with colleagues

The app needs Python packages **and** external binaries (`modkit`,
`samtools`/`htslib`). Two supported ways to distribute it:

### Option A — conda environment (recommended, incl. Apple Silicon)

Colleagues install [Miniforge](https://github.com/conda-forge/miniforge) once,
then, from a clone of this repo:

```bash
bash setup.sh                 # creates the "methyl-trio-ui" env + installs everything
conda activate methyl-trio-ui
methyl-trio                   # launches the app at http://localhost:8501
```

`setup.sh` builds the env from `environment.yml` (which pulls `ont-modkit`,
`samtools`, `htslib`, and — via `pip install -e .` — Streamlit, pysam, pandas,
scipy, matplotlib, plotly, and paramiko) and registers the `methyl-trio`
command. Equivalent manual steps:

```bash
conda env create -n methyl-trio-ui -f environment.yml
conda activate methyl-trio-ui
pip install -e .
methyl-trio
```

### Option B — Docker (zero local setup)

Colleagues with Docker need nothing else; the image bundles all tools:

```bash
export METHYL_TRIO_DATA=/absolute/path/to/data
docker compose up --build
```

### What each colleague does

1. Get the code: `git clone <repo-url>` (and `git checkout <branch>` if this is
   still on a feature branch).
2. Install with Option A or B above.
3. Launch and open <http://localhost:8501>.
4. In **Setup**, choose hg38/hg19 and click **Download and prepare** once (a
   one-time internet download, cached afterwards for offline use).
5. Select local BAMs (or connect to BCM on the **Remote data** page — each
   person logs in with their **own** credentials, which are never stored).

No shared secrets are required: references download per-machine into
`METHYL_TRIO_REFERENCE_CACHE`, and remote access uses each user's own login.

## Local inputs and references

- **Local file browser.** Setup selects each modBAM and the optional phased VCF
  from a rooted browser anchored at `METHYL_TRIO_DATA_ROOT` (`/data` under
  Docker, the user home directory natively). Click through folders, **paste a
  path into "Go to path"** to jump straight to a deep folder or file (server
  paths under a mount are translated automatically), and use the **filter** box
  to narrow large folders. Navigation cannot escape a data root.
- **External drives.** `METHYL_TRIO_DATA_ROOT` accepts several locations
  separated by the OS path separator (`:` on macOS/Linux, `;` on Windows), and a
  "Location" selector switches between them. When the variable is unset, the
  home directory plus any mounted external drives (`/Volumes`, `/media`,
  `/run/media`, `/mnt`) are browsable automatically. Under Docker, mount the
  drive read-only (e.g. `-v /Volumes/MyDrive:/data-external:ro`) and set
  `METHYL_TRIO_DATA_ROOT=/data:/data-external`.
- **BAM index detection.** The app looks for `<file>.bam.bai` or `<file>.bai`
  (also `.csi`) beside each BAM and explains how to create a missing index with
  `samtools index`.
- **Remote data over SSH (optional).** The "Remote data" page offers two ways
  to reach BAMs on a remote server, both using each person's own credentials
  held only in memory (never written to disk) over a two-hop `login1` →
  `analysis1` connection:
  - **Connect to BCM (in-app SSH):** log in with username/password (via
    `paramiko`), browse the server tree by path, and download the specific files
    you need into a local folder that then becomes browsable in Setup. No FUSE
    install required. `login1`'s host key is verified against `~/.ssh/known_hosts`.
    A **region-slice** action streams only the CpG/gene regions you name from a
    remote BAM (read-only `samtools view` on the server) into a small local,
    indexed BAM — no write permission on the server and only megabytes
    transferred, ideal for a targeted gene panel.
  - **Mount with SSHFS:** the app prints the exact `sshfs` command (with
    `ProxyJump`) for you to run in your own terminal — so your password/2FA stay
    between you and the server — then registers the read-only mount so BAMs are
    analyzed in place without a full download. Targeted runs are efficient over
    the mount; genome-wide runs are better run on the server. Requires SSHFS
    (macOS: macFUSE or FUSE-T + `sshfs`; Linux: `sshfs`).
- **Managed reference assemblies.** Pick hg38 or hg19; the FASTA, GENCODE GTF,
  and UCSC CpG islands are streamed once, written atomically, decompressed,
  FASTA-indexed with pysam, and cached under `METHYL_TRIO_REFERENCE_CACHE`
  (defaults to `~/.cache/methyl-trio/references`). Interrupted downloads are
  cleaned up and later runs reuse the cache offline.

## Workflow and method

1. Setup validates exactly three unique samples, one proband, BAM/reference
   contig lengths, reported basecaller models, and HP-tag availability. Optional
   relationship, clinical status, and phased-VCF fields include inline
   explanations and improve evidence-quality reporting.
2. Regions selects the analysis mode: genome-wide, selected chromosomes, or a
   GENCODE-derived promoter/gene-body gene panel.
3. Thresholds renders controls and rationale from `core/thresholds.py`.
4. Run executes the analysis and lets you choose where results are saved
   (the project `./runs` folder or any browsable data root, including an
   external drive).
5. Results provides a table, verdict, caveats, self-contained HTML, a one-click
   complete-run ZIP archive, and visualizations that mirror the reference
   scripts: per-gene % CpG methylation + coverage plots with the promoter
   highlighted and promoter/body heatmaps for targeted panels (script 07); and
   for genome-wide / per-chromosome runs a threshold sweep (each individual's
   private-DMR count vs effect size, with the sibling as the symmetric control),
   a karyotype-style DMR map with gene/imprinted labels, and an effect histogram
   against the family null (scripts 05/11). All figures are written to
   `<run>/figures/` and embedded in `report.html`. Results also includes an
   **interactive, zoomable** sliding-window methylation profile (Plotly): pick a
   gene/DMR or type `chrom:start-end`, drag to zoom, and pan with the range-slider
   overview track.

## Analysis methodology

The pipeline mirrors the reference trio scripts.

- **Pileup.** `modkit pileup --cpg --combine-strands --modified-bases 5mC
  --filter-threshold 0.7` per chromosome, then confident 5hmC calls (the
  `N_other` column) are folded into `N_canonical` so `valid_coverage ==
  N_mod + N_canonical` and `modkit dmr pair` accepts every row. Output is
  sorted, bgzipped, and tabix-indexed. (`--combine-strands` is toggleable in
  Setup for modBAMs without `MN` tags; the modified base is selectable.)
- **Genome-wide / selected chromosomes.** Per-chromosome `modkit dmr pair
  --segment` (HMM segmentation) against a one-chromosome FASTA to bound memory.
  A region is proband-private when it differs from **both** relatives in the
  same direction, with `num_sites` above the minimum and a Cohen's h 95% CI that
  excludes zero. The proband-vs-relative-2 / relative-1-vs-relative-2
  (relative–relative) comparison is the empirical null: its `|effect|`
  percentile sets the effect cutoff, and any region also variable there is
  subtracted. Known imprinted gDMRs are flagged as positive controls. chrX is
  only kept when both comparisons are between same-sex-female samples.
- **Targeted gene panel.** Promoter (TSS ± 2 kb) and gene-body (± 5 kb) regions
  are built from GENCODE, the panel is piled up, and each region is scored by a
  paired Wilcoxon signed-rank test on the CpGs covered at ≥ min-coverage in all
  three samples. A region is a candidate when it differs from both relatives by
  ≥ the minimum delta (percentage points) with `p < alpha`, concordant in
  direction, while the two relatives do not differ there. chrX regions are
  compared to the female relative only.

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

`pip install` provides the Python packages, but the pipeline also shells out to
the external binaries `modkit` and `tabix` (from `samtools`/`htslib`). Install
those from conda/bioconda — either create the full environment from
`environment.yml`, or add the tools to an existing environment. Note that
bioconda publishes modkit as `ont-modkit` (this also has an Apple Silicon
`osx-arm64` build):

```bash
conda env create -f environment.yml          # creates the "methyl-trio" env, or:
conda install -c bioconda -c conda-forge ont-modkit samtools htslib
```

```bash
python -m pip install -e '.[test]'
pytest
# Optional: point the browser and cache somewhere convenient.
export METHYL_TRIO_DATA_ROOT="$HOME/methylation-data"
export METHYL_TRIO_REFERENCE_CACHE="$HOME/.cache/methyl-trio/references"

# Launch (any one of these):
methyl-trio                 # console command installed by `pip install`
./run.sh                    # helper script in the repo root
streamlit run app/streamlit_app.py
```

`METHYL_TRIO_DATA_ROOT` bounds the Setup file browser; leave it unset to browse
from your home directory and any mounted external drives, or set it to one or
more `:`-separated locations to pin the browsable roots (for example
`export METHYL_TRIO_DATA_ROOT="$HOME/methylation-data:/Volumes/MyDrive"`).
`METHYL_TRIO_REFERENCE_CACHE` chooses where downloaded assemblies are stored so
they persist across runs.

## Outputs

Each run writes:

- `config.json` — resolved run manifest (including the selected assembly);
- `pipeline.log` — tool output;
- three pairwise DMR tables and sample bedMethyl files;
- `proband_specific_DMRs.tsv`;
- `dmr_effects.png`, `summary.json`, and `report.html`.

The Results page can bundle the entire run directory into a deterministic
`complete_run.zip` (manifests, logs, pileups, pairwise outputs, tables, figures,
and the HTML report, excluding the archive itself) with a local path and a
browser download button.

## Scientific limitations

This is a family prioritization screen, not a cohort analysis or diagnostic
test. Three related samples cannot estimate population variance. Basecalling,
sequencing, and processing batches can mimic methylation effects. Missing
parental lineages leave parent-of-origin mQTL effects unresolved; the report
generates a composition-specific warning. HP tags are optional for DMR calling
but required for haplotype read-level interpretation.

A phased VCF and identified mother and father make parent-of-origin and mQTL
follow-up possible, but their presence alone does not establish either effect.

## Validation status

Unit and synthetic smoke tests cover trio validation, sex-chromosome rules,
null-threshold intersection/ranking, `N_other` correction, command construction,
threshold bounds, targeted statistics, HTML export, safe data-root traversal,
BAM-index detection, reference manifests/cache reuse, interrupted-download
cleanup, tissue/batch metadata removal, complete-ZIP contents, and Streamlit
Setup/Results smoke tests.

BH11998 parity is not claimed by this repository: the original `00`–`08`
scripts, external-drive BAMs, and expected `proband_specific_DMRs.tsv` were not
present in the source repository. To perform the requested parity check, mount
those inputs, run the pilot chromosomes (`chr1`, `chr2`, `chr11`, `chr15`), and
compare coordinates/effects against the historical output.
