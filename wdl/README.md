# MANGO trio DMR workflow (WDL) for AnVIL / Terra

This is the **batch** route for running MANGO on AnVIL: a WDL workflow that runs
the exact same pipeline as the app via the `mango-run` CLI, on Cromwell, reading
`gs://` modBAMs that Terra localizes automatically. It is the better fit for
genome-wide runs on very large BAMs and for sharing a reproducible pipeline with
the GREGoR consortium through Dockstore. (For interactive, plot-driven analysis
of a single trio, use the MANGO app route — the interactive AnVIL Cloud
Environment integration — instead.)

## What it produces

`out/proband_specific_DMRs.tsv` (+`.bed` when candidates exist), per-comparison
segment/targeted tables, `summary.json`, `report.html`, figures, and a
`mango_run.tar.gz` of the whole run directory.

## Files

- [`mango_trio.wdl`](mango_trio.wdl) — the workflow + task.
- [`inputs.template.json`](inputs.template.json) — fill in `gs://` paths + the docker image.
- [`../.dockstore.yml`](../.dockstore.yml) — registers this workflow on Dockstore.

## Inputs at a glance

- Three samples (proband + two relatives): each takes a `*_bam`, `*_bai`, `*_sex`
  (`F`/`M`), and optional `*_affection` / `*_relationship`.
- `reference_fasta` **and** `reference_fai` (matching `.fai`), plus `gtf` (required
  for `mode = targeted`) and optional `cpg_islands`.
- `mode`: `targeted` (needs `genes`), `chromosomes` (needs `chromosomes`), or
  `whole_genome`.
- `docker`: your published MANGO image (see below).

---

# What YOU need to do (the manual, account-specific steps)

These steps require your credentials/accounts and cannot be done from the repo:

### 1. Build and publish the Docker image
The workflow needs a container that has `mango-run` plus `modkit`/`samtools`.
[`wdl/Dockerfile`](Dockerfile) builds exactly that (a minimal image; the app's
optional viz extras are excluded so it is small and builds reliably).

**Easiest (no cloud billing) — GitHub Actions -> GHCR.** This repo includes
`.github/workflows/build-image.yml`, which builds the image and pushes it to
GitHub Container Registry using the free built-in token. It runs automatically
when the branch is pushed, or you can trigger it from the repo's **Actions** tab
("Build and push MANGO image (GHCR)" -> **Run workflow**). After the first run,
make the package public so Terra can pull it anonymously: GitHub -> your profile
-> **Packages** -> `mango` -> **Package settings** -> **Change visibility** ->
**Public**. Then set:

```json
"mango_trio.docker": "ghcr.io/<your-github-username>/mango:latest"
```

**Alternative — Google Artifact Registry** (needs a GCP project with billing):

```bash
IMG=<REGION>-docker.pkg.dev/<YOUR_GCP_PROJECT>/<REPO>/mango:latest
docker build -f wdl/Dockerfile -t "$IMG" . && docker push "$IMG"
```

Then set `"mango_trio.docker": "<that image>"`. For a private Artifact Registry,
make sure your Terra proxy group has read access; a public GHCR/Quay image needs
no extra access.

### 2. Stage a reference bundle in your workspace bucket
Upload (once) an indexed FASTA plus GENCODE GTF and CpG islands to your bucket,
or point at an existing AnVIL reference bundle:

```bash
gsutil cp hg38.fa hg38.fa.fai gencode.v49.annotation.gtf cpgIslandExt.txt \
  gs://<your-workspace-bucket>/references/
```

Update the `reference_*`, `gtf`, and `cpg_islands` inputs to those `gs://` paths.

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
