#!/usr/bin/env nextflow

include { validateParameters; paramsHelp; paramsSummaryLog; fromSamplesheet } from 'plugin/nf-validation'


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
  exit 0
}
validateParameters()
log.info paramsSummaryLog(workflow)

file(params.out_dir + "/nextflow").mkdirs()
// groovy fonction within nextflow script
def parse_sample_entry(it) {
  def files = []
  if (it[1] &&!it[1].isEmpty()) files << file(it[1], checkIfExists: true)
  if (it[2] &&!it[2].isEmpty()) files << file(it[2], checkIfExists: true)
  def genome_ref = (it[3] &&!it[3].isEmpty())? file(it[3], checkIfExists: true) : null
  def kraken2_db = []
  if (it[4] &&!it[4].isEmpty()) {
    it[4].split(/;/).eachWithIndex() { db_path, idx -> 
      kraken2_db << file(db_path, checkIfExists: true)
    }
  }
  else {
    kraken2_db = null
  }
  meta = [
    "id": it[0],
    "assembly_type": it[5],
    "assembler": "spades", // make sure it's possible
    "realign": it[6]
  ]
  return [meta, files, kraken2_db, genome_ref]
}

def make_spades_args(meta) {
  switch (meta.assembly_type) {
    case "metaviral":
      return "--metaviral"
    case "corona":
      return "--corona"
    case "rnaviral":
      return "--rnaviral"
    default:
      if (meta.assembler == "spades") {
        error "Unknown value of 'assembly_type': ${meta.assembly_type}. Supported values are'metaviral', 'corona', and 'rnaviral'."
      } else {
        return null
      }
  }
}
  
// include
include {PRIMARY_FROM_READS} from './modules/subworkflows/primary/from_reads/main.nf'
include {VIRAL_ASSEMBLY} from './modules/subworkflows/assembly/viral_assembly/main.nf'

workflow {
  // START PARSING SAMPLE SHEET
  Channel.fromSamplesheet("input")
  | map {
    return parse_sample_entry(it)
  }
  | set { rawInputs }
  
  rawInputs
  | map { [it[0], it[1]] }
  | set {readsInputs}

  rawInputs
  | filter { it[2] }
  | map { [it[0].id, it[2]] }
  | set {k2Inputs}

  rawInputs
  | filter { it[3] }
  | map { [it[0].id, it[3]] }
  | set {refGenomeInputs}
  
  // START PRIMARY
  if (params.skip_primary) { // Note: skip_primary, include skiping spring_decompress step !!
    readsInputs
    | map {
      if (it[1].size() == 1) {
        it[0].read_type = "SR"
      } else {
        it[0].read_type = "PE"
      }
      return it
    }
    | set {trimmedInputs}
  } else {
    if ( params.kraken2_db == null ) {
      error "kraken2_db argument required for primary analysis"
    }

    Channel.fromPath(params.kraken2_db, type: "dir", checkIfExists: true)
    | collect
    | set {dbPathKraken2}
    
    PRIMARY_FROM_READS(readsInputs, dbPathKraken2)
    PRIMARY_FROM_READS.out.trimmed
    | set { trimmedInputs }
  }
  // END PRIMARY

  //// reads channel format: update meta
  trimmedInputs
  | map {
    it[0].args_spades = make_spades_args(it[0])
    return [it[0].id, it]
  }
  | join(k2Inputs, by: 0, remainder: true)
  | join(refGenomeInputs, by: 0, remainder: true)
  | set {inputsForViralAssembly}

  VIRAL_ASSEMBLY(inputsForViralAssembly)
  
}
