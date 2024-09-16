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
  log.info showSchemaHelp("assets/k2_db_schema.json")
  log.info showSchemaHelp("assets/ref_genome_schema.json")
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
  meta = [
    "id": it[0],
    "read_type": type,
    "ref_id": (it[3] && !it[3].isEmpty() ) ? it[3] : null,
    "k2_ids": (it[4] && !it[4].isEmpty() ) ? it[4].split(/;/) : [],
    "assembler": it[5].split(/;/),
    "realign": it[6]
  ]
  return [meta, files]
}
  
// include
include {PRIMARY_FROM_READS} from './modules/subworkflows/primary/from_reads/main.nf'
include {VIRAL_ASSEMBLY} from './modules/subworkflows/viral_assembly/main.nf'

workflow {
  // START PARSING SAMPLE SHEET
  Channel.fromSamplesheet("input")
  | map {
    return parse_sample_entry(it)
  }
  | set { readsInputs }

  Channel.fromSamplesheet("k2_dbs")
  | map { [["id": it[0]], it[1]] }
  | set {k2Inputs}

  Channel.fromSamplesheet("ref_genomes")
  | map { [["id": it[0]], it[1]] }
  | set {refGenomeInputs}

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

    numReads = Channel.value(params.num_reads_sample_qc)
    
    PRIMARY_FROM_READS(readsInputs, dbPathKraken2, numReads)
    PRIMARY_FROM_READS.out.trimmed
    | set { trimmedInputs }
  }
  // END PRIMARY

  VIRAL_ASSEMBLY(trimmedInputs, k2Inputs, refGenomeInputs)

  publish:
  VIRAL_ASSEMBLY.out.quast_dir            >> 'quast'
  VIRAL_ASSEMBLY.out.all_scaffolds        >> 'all_scaffolds'
  VIRAL_ASSEMBLY.out.all_aln              >> 'all_aln'
  PRIMARY_FROM_READS.out.trimmed          >> 'fastp'
  PRIMARY_FROM_READS.out.fastqc_trim_html >> 'fastqc_trim'
  PRIMARY_FROM_READS.out.fastqc_raw_html  >> 'fastqc_raw'
  PRIMARY_FROM_READS.out.multiqc_html     >> 'primary_multiqc'
}

output {
  directory "${params.out_dir}"
  mode params.publish_dir_mode
  'fastp' {
    enabled params.save_fastp
  }
  'all_aln' {
    enabled params.save_aln
  }
}
