#!/usr/bin/env nextflow
include {samplesheetToList; validateParameters} from 'plugin/nf-schema'
include {PRIMARY} from './modules/subworkflows/primary/main.nf'
include {VIRAL_ASSEMBLY} from './modules/subworkflows/viral_assembly/main.nf'
nextflow.preview.output = true

// groovy fonction within nextflow script
def parse_sample_entry(it) {
  def type = "SR"
  def files = [file(it[1])]
  if (it[2] && !it[2].isEmpty() ) {
    files << file(it[2])
    type = "PE"
  }
  if (it[1].toString().toLowerCase().endsWith("sfq")) {
    type = "sfq"
  }
  def meta = [
    "id": it[0],
    "read_type": type,
    "ref_id": (it[3] && !it[3].isEmpty() ) ? it[3] : null,
    "anchor_id": (it[4] && !it[4].isEmpty() ) ? it[4] : null,
    "class_db_ids": (it[5] && !it[5].isEmpty() ) ? it[5].split(/;/) : [],
    "class_tool": it[6],
    "assembler": it[7].split(/;/),
    "realign": it[8],
    "do_abacas": it[9],
    "keep_before_abacas": it[10],
    "dedup": it[11],
    "keep_before_dedup": it[12],
    "proteome_id": ""
  ]
  if (it[13] != "") {
    meta.proteome_id = it[13]
    meta.keep_before_hannot = it[14]
    meta.filter_annot = "yes"
    meta.revcomp = "yes"
    meta.retain_only_annot = "yes"
  }

  return [meta, files]
}

workflow {
  main:
  validateParameters()
  Channel.fromList(samplesheetToList(params.input, "assets/input_schema.json"))
  | map { it -> 
    return parse_sample_entry(it)
  }
  | set { readsInputs }

  if (params.class_dbs) {
    Channel.fromList(samplesheetToList(params.class_dbs, "assets/class_dbs_schema.json"))
    | map { it ->  [["id": it[0]], it[1]] }
    | set {k2Inputs}
  } else {
    k2Inputs = Channel.empty()
  }

  
  if (params.ref_genomes) {
    Channel.fromList(samplesheetToList(params.ref_genomes, "assets/ref_genomes_schema.json"))
    | map { it ->  [["id": it[0]], it[1]] }
    | set {refGenomeInputs}
  } else {
    refGenomeInputs = Channel.empty()
  }

  if (params.anchors) {
    Channel.fromList(samplesheetToList(params.anchors, "assets/anchors_schema.json"))
    | map { it ->  [["id": it[0]], it[1]] }
    | set {anchorsInputs}
  } else {
    anchorsInputs = Channel.empty()
  }

  if (params.prot) {
    Channel.fromList(samplesheetToList(params.prot, "assets/prot_schema.json"))
    | map { it -> [[id: it[0], regex_prot_name: it[2]], file(it[1])]}
    | set {protFasta}
  } else {
    protFasta = Channel.empty()
  }

  // START PRIMARY
  if (params.skip_primary) {
    trimmedInputs = readsInputs
  } else {
    if ( params.kraken2_db == null ) {
      error "kraken2_db argument required for primary analysis"
    }

    Channel.fromPath(params.kraken2_db, type: "dir", checkIfExists: true)
    | map { it -> [["id": "kraken_db"], it]}
    | collect
    | set {dbPathKraken2}

    taxDir = Channel.fromPath(params.tax_dir, type: 'dir')

    numReads = Channel.value(params.num_reads_sample_qc)
    
    PRIMARY(readsInputs, dbPathKraken2, taxDir, numReads)
    PRIMARY.out.trimmed
    | set { trimmedInputs }
  }

  VIRAL_ASSEMBLY(trimmedInputs, k2Inputs, refGenomeInputs, anchorsInputs, protFasta)

  publish:
  all_aln = VIRAL_ASSEMBLY.out.all_aln
  ref_aln = VIRAL_ASSEMBLY.out.ref_mapped_bam
  all_scaffolds = VIRAL_ASSEMBLY.out.all_scaffolds
  anchored_reads = VIRAL_ASSEMBLY.out.anchored_reads
  class_report = PRIMARY.out.class_report
  cleaned_reads = VIRAL_ASSEMBLY.out.cleaned_reads
  fastqc_for_raw = PRIMARY.out.fastqc_raw_html
  fastqc_for_trimmed = PRIMARY.out.fastqc_trim_html
  hannot_filtered = VIRAL_ASSEMBLY.out.hannot_filtered
  hannot_raw = VIRAL_ASSEMBLY.out.hannot_raw
  kraken2_report = PRIMARY.out.kraken2_report
  multiqc = PRIMARY.out.multiqc_html
  post_abacas_scaffolds = VIRAL_ASSEMBLY.out.post_abacas_scaffolds
  post_hannot_scaffolds = VIRAL_ASSEMBLY.out.post_hannot_scaffolds
  pre_abacas_scaffolds = VIRAL_ASSEMBLY.out.pre_abacas_scaffolds
  quast_dir = VIRAL_ASSEMBLY.out.quast_dir
  trimmed_and_filtered = PRIMARY.out.trimmed
  unclassed_reads = VIRAL_ASSEMBLY.out.unclassed_reads
}

output {
  all_aln {
    enabled params.save_aln
    path "mapping/vs_assemblies"
  }
  ref_aln {
    enabled params.save_aln
    path "mapping/vs_ref"
  }
  all_scaffolds {
    path "assembly/scaffolds"
  }
  anchored_reads {
    enabled params.save_anchored
    path "filtering/anchored_reads"
  }
  class_report {
    path "primary/classification"
  }
  cleaned_reads {
    enabled params.save_clean
    path "filtering/cleaned_reads"
  }
  fastqc_for_raw {
    path "primary/qc/fastqc/raw"
  }
  fastqc_for_trimmed {
    path "primary/qc/fastqc/trimmed"
  }
  hannot_filtered {
    path "assembly/annotation/hannot_filtered"
  }
  hannot_raw {
    path "assembly/annotation/hannot_raw"
  }
  kraken2_report {
    path "primary/classification/"
  }
  multiqc {
    path "primary/"
  }
  post_abacas_scaffolds {
    path "assembly/post_abacas_scaffolds"
  }
  post_hannot_scaffolds {
    path "assembly/post_hannot_scaffolds"
  }
  pre_abacas_scaffolds {
    path "assembly/pre_abacas_scaffolds"
  }
  quast_dir {
    path "assembly/quast"
  }
  trimmed_and_filtered {
    enabled params.save_fastp
    path "primary/trimmed_reads"
  }
  unclassed_reads {
    enabled params.save_unclassed
    path "filtering/unclassed_reads"
  }
}
