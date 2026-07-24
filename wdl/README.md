# MANGO trio DMR workflow (WDL) for AnVIL / Terra

This is the **batch** route for running MANGO on AnVIL: a WDL workflow that runs
the exact same pipeline as the app via the `mango-run` CLI, on Cromwell, reading
`gs://` modBAMs that Terra localizes automatically. It is the better fit for
genome-wide runs on very large BAMs and for sharing a reproducible pipeline with
the GREGoR consortium through Dockstore. (For interactive, plot-driven analysis
of a single trio, use the MANGO app route — the interactive AnVIL Cloud
Environment integration — instead.)

---

## Guideline: how consortium users run ONT vs PacBio

### Which Dockstore version?

Use workflow version **`v0.1.15`** (or newer `v*` / `main`).

Older tags (e.g. `v0.1.9`) lack the platform fix and recent modkit fixes — do not
use them for new runs.

### How do I choose Nanopore vs PacBio?

**Today (recommended):** import the published workflow **`mango_trio`**, then set
the **`platform`** input:

| Your data | Set `mango_trio.platform` to |
|-----------|------------------------------|
| Oxford Nanopore (Dorado) modBAMs | `"ont"` |
| PacBio HiFi modBAMs | `"pacbio"` |

That forces the correct modkit pileup settings. Do **not** rely on auto-detect
unless you know your BAMs carry Dorado `MN` tags (PacBio usually does not).

**Later (when published):** Dockstore will also expose dedicated workflows
`mango_ont` (MANGO-ONT) and `mango_pacbio` (MANGO-PACBIO) with the platform
baked in. Until those appear in public Dockstore search, use `mango_trio` +
`platform` as above. (Repo maintainers: publish them under Dockstore →
*My Workflows* if they are still unpublished.)

### Step-by-step in Terra

1. **Import:** Workflows → *Find a Workflow* → *Dockstore.org* → search
   `mango_trio` (path:
   `github.com/shaybehesht/methylation-pipeline/mango_trio`) → import into your
   workspace.
2. **Version:** open the workflow → select **`v0.1.15`** (not an old tag).
3. **Run mode:** choose “Run workflow(s) with inputs defined by file paths”
   (simplest for one trio).
4. **Fill in your data:**

   | Input | What to set |
   |-------|-------------|
   | `proband_bam`, `proband_bai` | proband `gs://` modBAM + `.bai` |
   | `relative1_bam/bai`, `relative2_bam/bai` | relatives’ modBAM + `.bai` |
   | `proband_sex`, `relative1_sex`, `relative2_sex` | `"F"` or `"M"` |
   | **`platform`** | **`"ont"` or `"pacbio"`** (required for correct analysis) |
   | `mode` | `"targeted"` (+ `genes`), `"chromosomes"` (+ `chromosomes`), or `"whole_genome"` |
   | `genes` / `chromosomes` | e.g. `["MECP2","SNRPN"]` / `["chr14"]` |
   | `disk_gb` | **`500`+** for WGS long-read BAMs (this is **disk**, not RAM) |
   | `memory_gb` | **`32`** (raise only if the job OOMs; do **not** set to 500) |
   | `preemptible` | `0` to avoid preemption restarts |
   | `docker` | see below |

5. **Leave at defaults unless you have a reason to change them:** `assembly`
   (`hg38`), blank `reference_fasta` / `gtf` / `cpg_islands` (auto-download).
6. **Optional:** `*_label`, `*_affection`, `*_relationship` improve the report.
7. **Launch.** Outputs land under the submission’s `call-run_mango/` folder
   (`report.html`, `summary.json`, tables, figures, `mango_run.tar.gz`).

### Docker image

Always include the `docker.io/` prefix (Terra/Cromwell requirement):

```text
docker.io/shaghayeghb/mango:latest
```

To pin a known-good build (recommended for papers / shared runs):

```text
docker.io/shaghayeghb/mango:870245c16f21be0c711deae379e7270e4a075cff
```

(`870245c…` is the v0.1.14+ image that includes ONT/PACBIO platform support and
modkit fixes. Prefer `latest` only if you are okay tracking new builds.)

### Memory vs disk (common mix-up)

| Input | Typical value | Meaning |
|-------|----------------|---------|
| `disk_gb` | **500** | Scratch disk for BAMs + reference download |
| `memory_gb` | **32** | RAM for `mango-run` / modkit |

### Requirements / caveats

- You need your own Terra billing project; runs bill to you.
- Inputs must be **modBAM + `.bai`** (CRAM is not supported — convert first).
- Do not mix ONT and PacBio within one trio.
- Empty `proband_specific_DMRs.tsv` with a successful run can mean “no private
  candidates passed filters,” not a crash — check `summary.json`
  (`candidate_count` vs `null_count`) and `*.segments.tsv`.

### What “version” means

The Dockstore/Terra **version dropdown** (e.g. `v0.1.15`, `main`) is the WDL
snapshot. **ONT vs PacBio is not a version** — it is either:

- the `platform` input on `mango_trio`, or
- a separate workflow (`mango_ont` / `mango_pacbio`) once those are published.

---

## What it produces

`out/proband_specific_DMRs.tsv` (+`.bed` when candidates exist), per-comparison
segment/targeted tables, `summary.json`, `report.html`, figures, and a
`mango_run.tar.gz` of the whole run directory.

## Files

- [`mango_ont.wdl`](mango_ont.wdl) — **MANGO-ONT** workflow (Dorado / Nanopore).
- [`mango_pacbio.wdl`](mango_pacbio.wdl) — **MANGO-PACBIO** workflow (PacBio HiFi).
- [`mango_trio.wdl`](mango_trio.wdl) — published workflow; set `platform` to `ont` or `pacbio`.
- [`mango_run_task.wdl`](mango_run_task.wdl) — shared Cromwell task.
- [`inputs.ont.template.json`](inputs.ont.template.json) / [`inputs.pacbio.template.json`](inputs.pacbio.template.json) / [`inputs.template.json`](inputs.template.json) — input templates.
- [`../.dockstore.yml`](../.dockstore.yml) — registers workflows on Dockstore (`main` + `v*` tags only).

## Inputs at a glance

- Three samples (proband + two relatives): each takes a `*_bam`, `*_bai`, `*_sex`
  (`F`/`M`), and optional `*_affection` / `*_relationship`.
- **`platform`:** `ont` | `pacbio` | `auto` (prefer `ont` or `pacbio`).
- Reference: either set `assembly` (`hg38`/`hg19`) to auto-download, or provide
  `reference_fasta` + `reference_fai` (+ `gtf` for `mode = targeted`, optional
  `cpg_islands`) to override.
- `mode`: `targeted` (needs `genes`), `chromosomes` (needs `chromosomes`), or
  `whole_genome`.
- `docker`: public MANGO image (see above).
- `disk_gb` / `memory_gb`: see table above.

---

# What YOU need to do (the manual, account-specific steps)

These steps require your credentials/accounts and cannot be done from the repo:

### 1. Build and publish the Docker image
The workflow needs a container that has `mango-run` plus `modkit`/`samtools`.
[`wdl/Dockerfile`](Dockerfile) builds exactly that (a minimal image; the app's
optional viz extras are excluded so it is small and builds reliably).

**Important:** Terra/Cromwell does **not** support `ghcr.io` for WDL execution
(the job fails with "Registry ghcr.io is not supported"). Use **Docker Hub**
(`docker.io`), which Terra supports and is free with no cloud billing.

**GitHub Actions -> Docker Hub.** `.github/workflows/build-image.yml` builds the
image and pushes it to Docker Hub when these two repo secrets are set
(GitHub repo -> Settings -> Secrets and variables -> Actions):

- `DOCKERHUB_USERNAME` - your Docker Hub username
- `DOCKERHUB_TOKEN` - a Docker Hub access token (Docker Hub -> Account Settings
  -> Personal access tokens -> Generate, read/write)

Create a free Docker Hub account (no credit card) at https://hub.docker.com,
add the secrets, then push to `main` (or use the **Actions** tab -> **Run
workflow**). The image is pushed to `docker.io/<username>/mango:latest`
(Docker Hub repos are public by default). Then set:

```json
"mango_trio.docker": "docker.io/<your-dockerhub-username>/mango:latest"
```

**Alternative — Google Artifact Registry** (needs a GCP project with billing):

```bash
IMG=<REGION>-docker.pkg.dev/<YOUR_GCP_PROJECT>/<REPO>/mango:latest
docker build -f wdl/Dockerfile -t "$IMG" . && docker push "$IMG"
```

Then set `"mango_trio.docker": "<that image>"`. For a private Artifact Registry,
make sure your Terra proxy group has read access; a public GHCR/Quay image needs
no extra access.

### 2. Reference — auto-download (default) or bring your own
**Easiest:** leave `reference_fasta`/`reference_fai`/`gtf`/`cpg_islands` unset and
set `assembly` to `hg38` (or `hg19`). The task downloads and prepares the matching
FASTA, GENCODE GTF, and CpG islands automatically (needs network egress from the
task, which Terra allows by default). No staging required.

**Bring your own (fully hermetic):** stage an indexed FASTA plus GENCODE GTF and
CpG islands once, then set the inputs to those `gs://` paths (they override the
download):

```bash
gsutil cp hg38.fa hg38.fa.fai gencode.v49.annotation.gtf cpgIslandExt.txt \
  gs://<your-workspace-bucket>/references/
```

Then set `reference_fasta`, `reference_fai`, `gtf`, and `cpg_islands` accordingly.

### 3. Get your BAM `gs://` paths from the workspace
Find the modBAMs (and their `.bai` indexes) in your GREGoR/AnVIL Data Table or
bucket, and put them in the inputs JSON (or bind to Data Table columns like
`this.proband_bam` when launching from a table). Requester-pays reads are billed
to your workspace's billing project automatically when run in Terra.

### 4. Run it in Terra (either way)
- **Import to a workspace:** in Terra → Workflows → *Find a Workflow* →
  *Import from Dockstore* (after step 6), or *Add a Workflow* → upload
  `mango_trio.wdl` directly. Provide the inputs JSON and **Run**.
- **Local validation (optional):** `womtool validate wdl/mango_trio.wdl` or
  `miniwdl check wdl/mango_trio.wdl`.

### 5. (If sharing with the consortium) Register on Dockstore
- Connect your GitHub repo to Dockstore; it will pick up `.dockstore.yml`.
- Publish the `mango_trio` workflow.
- **Contact the GREGoR DCC** to add the workflow/repo to the **GREGoR
  organization** on Dockstore so consortium members can import it. (This
  org-add step can only be done by the DCC.)

### 6. Verify a small run first
Run one trio in `targeted` mode on a few genes (fast, cheap) before launching
`whole_genome`, and confirm `report.html` / `proband_specific_DMRs.tsv` look
right.

---

## What is already done for you in this repo
- `mango-run` headless CLI (`core/cli.py`) and its console-script entry point.
- The WDL, inputs template, and `.dockstore.yml`.
- The Docker image definition (`wdl/Dockerfile`) that includes `mango-run`.
- A GitHub Actions workflow (`.github/workflows/build-image.yml`) that builds and
  pushes that image to GHCR with no cloud billing.
