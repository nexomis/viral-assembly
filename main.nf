#!/usr/bin/env nextflow

include { validateParameters; paramsHelp; paramsSummaryLog } from 'plugin/nf-validation'

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
  log.info paramsHelp("nextflow run nexomis/viral-assembly [args]")
  exit 0
}
validateParameters()
log.info paramsSummaryLog(workflow)

workflow {

/*
Define workflow here
To access to a parameter (default or mandatory) call param.name_param
*/


}