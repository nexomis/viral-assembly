#!/usr/bin/env nextflow
nextflow.preview.output = true
include { validateParameters; paramsHelp; paramsSummaryLog; fromSamplesheet } from 'plugin/nf-validation'
include { showSchemaHelp; extractType } from './modules/config/schema_helper.nf'

log.info """
    |            #################################################
    |            #    _  _                             _         #
    |            #   | \\| |  ___  __ __  ___   _ __   (_)  __    #
    |            #   | .` | / -_) \\ \\ / / _ \\ | '  \\  | | (_-<   #
    |            #   |_|\\_| \\___| /_\\_\\ \\___/ |_|_|_| |_| /__/   #
    |            #                                               #
    |            #################################################
    |
    | viral-assembly: Assemble viral genome and perform the associated QC. Starting from raw reads and one pre-built reference database in standard format (cf. README.md).
    |
    |""".stripMargin()

if (params.help) {
  log.info paramsHelp("nextflow run nexomis/viral-assembly --input </path/to/samplesheet> [args]")
  log.info showSchemaHelp("assets/input_schema.json")
  log.info showSchemaHelp("assets/class_dbs_schema.json")
  log.info showSchemaHelp("assets/ref_genomes_schema.json")
  exit 0
}
validateParameters()
log.info paramsSummaryLog(workflow)

file(params.out_dir + "/nextflow").mkdirs()
// groovy fonction within nextflow script
def parse_sample_entry(it) {
  def type = "SR"
  def files = [file(it[1])]
  if (it[2] && !it[2].isEmpty() ) {
    files << file(it[2])
    type = "PE"
  } else {
    if (it[1].toString().toLowerCase().endsWith("spring")) {
      type = "spring"
    }
  }
  def meta = [
    "id": it[0],
    "read_type": type,
    "ref_id": (it[3] && !it[3].isEmpty() ) ? it[3] : null,
    "class_db_ids": (it[4] && !it[4].isEmpty() ) ? it[4].split(/;/) : [],
    "class_tool": it[5],
    "assembler": it[6].split(/;/),
    "realign": it[7],
    "do_abacas": it[8],
    "keep_before_abacas": it[9],
    "dedup": it[10],
    "keep_before_dedup": it[11],
    "proteome_id": ""
  ]
  if (it[12] != "") {
    meta.proteome_id = it[12]
    meta.keep_before_hannot = it[13]
    meta.filter_annot = "yes"
    meta.revcomp = "yes"
    meta.retain_only_annot = "yes"
  }

  return [meta, files]
}
  
// include
include {PRIMARY} from './modules/subworkflows/primary/main.nf'
include {VIRAL_ASSEMBLY} from './modules/subworkflows/viral_assembly/main.nf'

workflow {
  // START PARSING SAMPLE SHEET
  Channel.fromSamplesheet("input")
  | map {
    return parse_sample_entry(it)
  }
  | set { readsInputs }

  if (params.class_dbs) {
    Channel.fromSamplesheet("class_dbs")
    | map { [["id": it[0]], it[1]] }
    | set {k2Inputs}
  } else {
    k2Inputs = Channel.empty()
  }

  
  if (params.ref_genomes) {
    Channel.fromSamplesheet("ref_genomes")
    | map { [["id": it[0]], it[1]] }
    | set {refGenomeInputs}
  } else {
    refGenomeInputs = Channel.empty()
  }

  if (params.prot) {
    Channel.fromSamplesheet("prot")
    | map {[[id: it[0], regex_prot_name: it[2]], file(it[1])]}
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
    | map {[["id": "kraken_db"], it]}
    | collect
    | set {dbPathKraken2}

    taxDir = Channel.fromPath(params.tax_dir, type: 'dir')

    numReads = Channel.value(params.num_reads_sample_qc)
    
    PRIMARY(readsInputs, dbPathKraken2, taxDir, numReads)
    PRIMARY.out.trimmed
    | set { trimmedInputs }
  }
  // END PRIMARY

  VIRAL_ASSEMBLY(trimmedInputs, k2Inputs, refGenomeInputs, protFasta)

  publish:
  PRIMARY.out.trimmed                         >> 'trimmed_and_filtered'
  PRIMARY.out.fastqc_trim_html                >> 'fastqc_for_trimmed'
  PRIMARY.out.fastqc_raw_html                 >> 'fastqc_for_raw'
  PRIMARY.out.multiqc_html                    >> 'multiqc'
  PRIMARY.out.kraken2_report                  >> 'classification'
  PRIMARY.out.class_report                    >> 'classification'
  VIRAL_ASSEMBLY.out.cleaned_reads            >> 'cleaned_reads'
  VIRAL_ASSEMBLY.out.quast_dir                >> 'quast'
  VIRAL_ASSEMBLY.out.all_scaffolds            >> 'all_scaffolds'
  VIRAL_ASSEMBLY.out.all_aln                  >> 'all_aln'
  VIRAL_ASSEMBLY.out.pre_abacas_scaffolds     >> 'pre_abacas'
  VIRAL_ASSEMBLY.out.post_abacas_scaffolds    >> 'post_abacas'
  VIRAL_ASSEMBLY.out.post_hannot_scaffolds    >> 'post_hannot'
  VIRAL_ASSEMBLY.out.hannot_raw               >> 'hannot_raw'
  VIRAL_ASSEMBLY.out.hannot_filtered          >> 'hannot_filtered'
}

output {
  directory "${params.out_dir}"
  mode params.publish_dir_mode
  'trimmed_and_filtered' {
    enabled params.save_fastp
  }
  'cleaned_reads' {
    enabled params.save_clean
  }
  'all_aln' {
    enabled params.save_aln
  }
}
