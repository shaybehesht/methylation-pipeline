version 1.0

# MANGO-PACBIO — trio methylation DMR for PacBio HiFi modBAMs.
# Forces --platform pacbio: general workers (--motif CG 0), no --combine-strands.

import "mango_run_task.wdl" as mango

workflow mango_pacbio {
  input {
    File proband_bam
    File proband_bai
    String proband_label = "proband"
    String proband_sex
    String? proband_affection

    File relative1_bam
    File relative1_bai
    String relative1_label = "relative1"
    String relative1_sex
    String? relative1_relationship
    String? relative1_affection

    File relative2_bam
    File relative2_bai
    String relative2_label = "relative2"
    String relative2_sex
    String? relative2_relationship
    String? relative2_affection

    File? reference_fasta
    File? reference_fai
    File? gtf
    File? cpg_islands
    String assembly = "hg38"

    String mode = "targeted"
    Array[String] chromosomes = []
    Array[String] genes = []

    File? phased_vcf
    Array[String] modified_bases = ["5mC"]
    Array[String] threshold_overrides = []

    String docker = "docker.io/shaghayeghb/mango:latest"
    Int cpu = 4
    Int memory_gb = 32
    Int disk_gb = 100
    Int preemptible = 1
  }

  call mango.run_mango {
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
      platform = "pacbio", combine_strands = false,
      threshold_overrides = threshold_overrides,
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
    description: "MANGO-PACBIO: trio methylation DMR for PacBio HiFi modBAMs (AnVIL/Terra)."
    author: "MANGO"
  }
}
