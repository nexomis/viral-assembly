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

def parse_sample_entry_wReadType(it) {  // sample_name,path_raw_R1,path_raw_R2,ref_genome,host_taxid,assembly_type(meta/viral/coronavirus)
// replaced by 'parse_sample_entry()' because 'read_type' is defined here and on SPRING_DECOMPRESS process (called by PRIMARY_FROM_READS worflow) ... 
  def meta = [ "id": it[0], "ref_genome": it[3], "host_txid": it[4], "assembly_type": it[5] ]
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

def parse_sample_entry(it) {  // sample_name,path_raw_R1,path_raw_R2,ref_genome,host_kraken_db,assembly_type(meta/viral/coronavirus)
  def meta = [ "id": it[0], "ref_genome": it[3], "host_kraken_db": it[4], "assembly_type": it[5] ]
  def r1_file = file(it[1])
  if (it[2] == "") {
    return [meta, [r1_file] ]
  } else {
    def r2_file = file(it[2])
    return [meta, [r1_file, r2_file] ]
  }
}

def make_spades_args(meta) {
  spades_args = []
  if (meta.assembly_type == "metaviral") {
    spades_args << "--metaviral"
  } else if (meta.assembly_type == "corona") {
    spades_args << "--corona"
  } else if (meta.assembly_type == "rnaviral") {
    spades_args << "--rnaviral"
  } else {
    spades_args << "meta.assembly_type"
    println "Unknown value of 'assembly_type': ${meta.assembly_type}. Argument integrated as is in the assembler command line!!!"
  }
  return spades_args.join(" ")
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
  // END PARSING SAMPLE SHEET

  // START PRIMARY
  if (params.skip_primary) { // Note: skip_primary, include skiping spring_decompress step !!
    rawInputs
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
    
    PRIMARY_FROM_READS(rawInputs, dbPathKraken2)
    
    trimmedInputs = PRIMARY_FROM_READS.out.trimmed
  }
  // END PRIMARY


  // START ASSEMBLY: KRAKEN2_SPADES_ABACAS_BOWTIE2wBUILD_QUAST

  //// reads channel format: update meta
  // TODO: match between 'meta.kraken_db_id of' and 'meta.id' of krakenDb ?!
  trimmedInputs
  | map {
    it[0].spades_args = make_spades_args(it[0])
    return it
  }
  | set {trimmedInputsWithAssemblyArgs}
  
  VIRAL_ASSEMBLY(trimmedInputsWithAssemblyArgs) 
  // END ASSEMBLY: KRAKEN2_SPADES_ABACAS_BOWTIE2wBUILD_QUAST

}
