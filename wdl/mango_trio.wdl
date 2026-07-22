version 1.0

# MANGO — three-sample (trio) nanopore methylation DMR analysis as a batch
# workflow for AnVIL / Terra (Cromwell). It runs the exact same pipeline as the
# MANGO app via the `mango-run` headless CLI, reading gs:// modBAMs that Cromwell
# localizes automatically and writing DMR tables, an HTML report, and figures.
#
# Share it with the GREGoR consortium by publishing this descriptor on Dockstore
# (see wdl/README.md).
#
# References: either provide reference_fasta (+ reference_fai, and gtf for
# targeted mode) explicitly for a fully hermetic run, OR leave them unset and
# set `assembly` (hg38/hg19) so the task downloads and prepares the FASTA, GTF,
# and CpG islands automatically (requires network egress from the task).

workflow mango_trio {
  input {
    # --- Proband ---
    File proband_bam
    File proband_bai
    String proband_label = "proband"
    String proband_sex               # "F" or "M"
    String? proband_affection        # affected | unaffected | unknown

    # --- Relative 1 ---
    File relative1_bam
    File relative1_bai
    String relative1_label = "relative1"
    String relative1_sex
    String? relative1_relationship   # mother | father | sibling | other
    String? relative1_affection

    # --- Relative 2 ---
    File relative2_bam
    File relative2_bai
    String relative2_label = "relative2"
    String relative2_sex
    String? relative2_relationship
    String? relative2_affection

    # --- Reference ---
    # Provide reference_fasta (+ reference_fai; gtf required for targeted mode),
    # OR leave all four unset and rely on `assembly` to auto-download them.
    File? reference_fasta
    File? reference_fai
    File? gtf
    File? cpg_islands
    String assembly = "hg38"          # used to auto-download when reference_fasta is unset

    # --- Region scope ---
    String mode = "targeted"         # whole_genome | chromosomes | targeted
    Array[String] chromosomes = []   # for mode = chromosomes
    Array[String] genes = []         # for mode = targeted

    # --- Options ---
    File? phased_vcf
    Array[String] modified_bases = ["5mC"]
    Boolean combine_strands = true
    Array[String] threshold_overrides = []   # e.g. ["min_sites=5", "alpha=0.05"]

    # --- Runtime ---
    String docker = "docker.io/shaghayeghb/mango:latest"
    Int cpu = 4
    Int memory_gb = 32
    Int disk_gb = 100
    Int preemptible = 1
  }

  call run_mango {
    input:
      proband_bam = proband_bam, proband_bai = proband_bai,
      proband_label = proband_label, proband_sex = proband_sex,
      proband_affection = proband_affection,
      relative1_bam = relative1_bam, relative1_bai = relative1_bai,
      relative1_label = relative1_label, relative1_sex = relative1_sex,
      relative1_relationship = relative1_relationship,
      relative1_affection = relative1_affection,
      relative2_bam = relative2_bam, relative2_bai = relative2_bai,
      relative2_label = relative2_label, relative2_sex = relative2_sex,
      relative2_relationship = relative2_relationship,
      relative2_affection = relative2_affection,
      reference_fasta = reference_fasta, reference_fai = reference_fai,
      gtf = gtf, cpg_islands = cpg_islands, assembly = assembly,
      mode = mode, chromosomes = chromosomes, genes = genes,
      phased_vcf = phased_vcf, modified_bases = modified_bases,
      combine_strands = combine_strands, threshold_overrides = threshold_overrides,
      docker = docker, cpu = cpu, memory_gb = memory_gb,
      disk_gb = disk_gb, preemptible = preemptible
  }

  output {
    File proband_specific_dmrs = run_mango.proband_specific_dmrs
    File? proband_specific_dmrs_bed = run_mango.proband_specific_dmrs_bed
    File summary = run_mango.summary
    File report = run_mango.report
    File run_config = run_mango.run_config
    File results_archive = run_mango.results_archive
    Array[File] tables = run_mango.tables
    Array[File] figures = run_mango.figures
  }

  meta {
    description: "Trio nanopore methylation DMR screen (MANGO) for AnVIL/Terra."
    author: "MANGO"
  }
}

task run_mango {
  input {
    File proband_bam
    File proband_bai
    String proband_label
    String proband_sex
    String? proband_affection

    File relative1_bam
    File relative1_bai
    String relative1_label
    String relative1_sex
    String? relative1_relationship
    String? relative1_affection

    File relative2_bam
    File relative2_bai
    String relative2_label
    String relative2_sex
    String? relative2_relationship
    String? relative2_affection

    File? reference_fasta
    File? reference_fai
    File? gtf
    File? cpg_islands
    String assembly

    String mode
    Array[String] chromosomes
    Array[String] genes

    File? phased_vcf
    Array[String] modified_bases
    Boolean combine_strands
    Array[String] threshold_overrides

    String docker
    Int cpu
    Int memory_gb
    Int disk_gb
    Int preemptible
  }

  command <<<
    set -euo pipefail
    mkdir -p refs bams out

    # Stage the BAMs with the sibling-index names the tools expect.
    ln -s "~{proband_bam}"     "bams/~{proband_label}.bam"
    ln -s "~{proband_bai}"     "bams/~{proband_label}.bam.bai"
    ln -s "~{relative1_bam}"   "bams/~{relative1_label}.bam"
    ln -s "~{relative1_bai}"   "bams/~{relative1_label}.bam.bai"
    ln -s "~{relative2_bam}"   "bams/~{relative2_label}.bam"
    ln -s "~{relative2_bai}"   "bams/~{relative2_label}.bam.bai"

    # Reference: either the provided FASTA (+.fai) or auto-download via --assembly.
    # Use default="" so the optional File placeholders render cleanly when unset.
    REFERENCE_FASTA="~{default="" reference_fasta}"
    REFERENCE_FAI="~{default="" reference_fai}"
    if [ -n "$REFERENCE_FASTA" ]; then
      ln -s "$REFERENCE_FASTA" refs/reference.fa
      ln -s "$REFERENCE_FAI"   refs/reference.fa.fai
      REF_ARG="--reference-fasta refs/reference.fa ~{"--gtf " + gtf} ~{"--cpg-islands " + cpg_islands}"
    else
      REF_ARG="--assembly ~{assembly} ~{"--gtf " + gtf} ~{"--cpg-islands " + cpg_islands}"
    fi

    # Build multi-value flags in bash from write_lines files. This avoids
    # Cromwell's "Cannot interpolate Array[Nothing]" error that sep= raises on an
    # empty array (e.g. genes=[] in chromosomes mode, or chromosomes=[] otherwise).
    CHROMS="$(tr '\n' ' ' < ~{write_lines(chromosomes)})"
    GENES="$(tr '\n' ' ' < ~{write_lines(genes)})"
    if [ -n "${CHROMS// /}" ]; then CHROM_ARG="--chromosomes $CHROMS"; else CHROM_ARG=""; fi
    if [ -n "${GENES// /}" ]; then GENE_ARG="--genes $GENES"; else GENE_ARG=""; fi

    MOD_ARGS=""
    while read -r _mb; do [ -n "$_mb" ] && MOD_ARGS="$MOD_ARGS --modified-base $_mb"; done < ~{write_lines(modified_bases)}
    THRESH_ARGS=""
    while read -r _th; do [ -n "$_th" ] && THRESH_ARGS="$THRESH_ARGS --set-threshold $_th"; done < ~{write_lines(threshold_overrides)}

    on_fail() {
      ec=$?
      echo "mango-run failed (exit ${ec})" >&2
      if [ -f out/pipeline.log ]; then
        echo "=== tail out/pipeline.log ===" >&2
        tail -n 80 out/pipeline.log >&2
      fi
      exit "${ec}"
    }
    trap on_fail ERR

    mango-run \
      --proband-bam "bams/~{proband_label}.bam" --proband-label "~{proband_label}" \
      --proband-sex ~{proband_sex} ~{"--proband-affection " + proband_affection} \
      --relative1-bam "bams/~{relative1_label}.bam" --relative1-label "~{relative1_label}" \
      --relative1-sex ~{relative1_sex} \
      ~{"--relative1-relationship " + relative1_relationship} \
      ~{"--relative1-affection " + relative1_affection} \
      --relative2-bam "bams/~{relative2_label}.bam" --relative2-label "~{relative2_label}" \
      --relative2-sex ~{relative2_sex} \
      ~{"--relative2-relationship " + relative2_relationship} \
      ~{"--relative2-affection " + relative2_affection} \
      $REF_ARG \
      --mode ~{mode} \
      $CHROM_ARG $GENE_ARG \
      ~{"--phased-vcf " + phased_vcf} \
      $MOD_ARGS \
      ~{if combine_strands then "" else "--no-combine-strands"} \
      $THRESH_ARGS \
      --output-dir out

    tar -czf mango_run.tar.gz -C out .
  >>>

  output {
    File proband_specific_dmrs = "out/proband_specific_DMRs.tsv"
    File? proband_specific_dmrs_bed = "out/proband_specific_DMRs.bed"
    File summary = "out/summary.json"
    File report = "out/report.html"
    File run_config = "out/config.json"
    File results_archive = "mango_run.tar.gz"
    Array[File] tables = glob("out/*.tsv")
    Array[File] figures = glob("out/figures/*.png")
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} SSD"
    preemptible: preemptible
  }
}
