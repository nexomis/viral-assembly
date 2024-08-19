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

def parse_sample_entry(it) {  // sample_name,path_raw_R1,path_raw_R2,ref_genome,host_taxid,type(meta/viral/coronavirus)
  def meta = [ "id": it[0], "ref_genome": it[3], "host_txid": it[4], "sample_type": it[5] ]
  def r1_file = file(it[1])
  if (it[2] == "") {
    if (r1_file.getExtension() == "spring") {
      meta.read_type = "spring"
    } else {
      meta.read_type = "SR"
    }
    return [meta, [r1_file] ]
  } else {
    def r2_file = file(it[2])
    meta.read_type = "PE"
    return [meta, [r1_file, r2_file] ]
  }
}

def make_spades_args(meta) {
  spades_args = []
  if (meta.sample_type == "metaviral") {
    spades_args << "--metaviral"
  } else if (meta.sample_type == "corona") {
    spades_args << "--corona"
  } else if (meta.sample_type == "rnaviral") {
    spades_args << "--rnaviral"
  } else {
    spades_args << "meta.sample_type"
    println "Unknown value of 'sample_type': ${meta.sample_type}. Argument integrated as is in the assembler command line!!!"
  }
  return spades_args.join(" ")
}

  
// include
include {SPRING_DECOMPRESS} from './modules/process/spring/decompress/main.nf'
include {PRIMARY_FROM_READS} from './modules/subworkflows/primary/from_reads/main.nf'
include {SPADES_ABACAS} from './modules/subworkflows/assembly/spades_abacas/main.nf'

workflow {
  // START PARSING SAMPLE SHEET
  Channel.fromSamplesheet("input")
  | map {
    return parse_sample_entry(it)
  }
  | branch {
    spring: it[0].read_type == "spring"
    fastq: it[0].read_type != "spring"
  }
  | set { inputs }

  SPRING_DECOMPRESS(inputs.spring)
  | map {
    if (it[1].size()==1) {
      it[0].read_type = "SR"
    } else {
      it[0].read_type = "PE"
    }
    return it
  }
  | concat(inputs.fastq)
  | set {rawInputs}
  // END PARSING SAMPLE SHEET

  // START PRIMARY
  if (params.skip_primary) {
    trimmedInputs = rawInputs
  } else {
    if ( params.kraken2_db == null ) {
      error "kraken2_db argument required for primary analysis"
    }
    Channel.fromPath(params.kraken2_db, type: "dir", checkIfExists: true)
    | collect
    | set {dbPathKraken2}
    
    PRIMARY_FROM_READS(rawInputs, dbPathKraken2)
    
    trimmedInputs = PRIMARY_FROM_READS.out.trimmed
  }
  // END PRIMARY

  // START ASSEMBLY:SPADES_ABACAS

  //// reads channel format: update meta
  trimmedInputs
  | map {
    it[0].spades_args = make_spades_args(it[0])
    return it
  }
  | set {inputsForSpadesAbacas}
  
  //// make abacas argument as channel
  if (params.abacas_MUMmer_program == "nucmer" || params.abacas_MUMmer_program == "promer") {
    abacas_MUMmer_program = params.abacas_MUMmer_program
  } else {
    abacas_MUMmer_program = "nucmer"
    println "Unknown value of 'abacas_MUMmer_program': '${params.abacas_MUMmer_program}'. Replaced by 'nucmer' value."
  }

  abacas_keep_on_output = params.abacas_keep_on_output ?: 'scaffold'

  abacasMUMmerProgram = Channel.value(abacas_MUMmer_program)
  abacasKeepOnOutput = Channel.value(abacas_keep_on_output)

  SPADES_ABACAS(inputsForSpadesAbacas, abacasMUMmerProgram, abacasKeepOnOutput)
  // END ASSEMBLY:SPADES_ABACAS

}
