version 1.0

# Shared Cromwell task for MANGO trio runs. Imported by mango_ont.wdl,
# mango_pacbio.wdl, and mango_trio.wdl (auto platform).

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
    String platform
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
    REFERENCE_FASTA="~{default="" reference_fasta}"
    REFERENCE_FAI="~{default="" reference_fai}"
    if [ -n "$REFERENCE_FASTA" ]; then
      ln -s "$REFERENCE_FASTA" refs/reference.fa
      ln -s "$REFERENCE_FAI"   refs/reference.fa.fai
      REF_ARG="--reference-fasta refs/reference.fa ~{"--gtf " + gtf} ~{"--cpg-islands " + cpg_islands}"
    else
      REF_ARG="--assembly ~{assembly} ~{"--gtf " + gtf} ~{"--cpg-islands " + cpg_islands}"
    fi

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
      if [ -f out/run.log ]; then
        echo "=== tail out/run.log ===" >&2
        tail -n 80 out/run.log >&2
      fi
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
      --platform ~{platform} \
      ~{if combine_strands then "" else "--no-combine-strands"} \
      $THRESH_ARGS \
      --output-dir out 2>&1 | tee out/run.log

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
