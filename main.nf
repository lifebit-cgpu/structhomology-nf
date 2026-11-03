#!/usr/bin/env nextflow
/*
 * structhomology-nf
 * Structure (Foldseek) versus sequence (MMseqs2) homology search for AlphaFold DB
 * models of UniProt accessions, over one shared target database.
 * https://github.com/lifebit-cgpu/structhomology-nf
 */

include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'
include { STRUCTHOMOLOGY } from './workflows/structhomology'

workflow {
    // Fails with a readable message if any parameter breaks nextflow_schema.json
    validateParameters()
    log.info paramsSummaryLog(workflow)

    // Samplesheet validated against assets/schema_input.json on entry:
    // a malformed accession, an unknown group or a duplicated row stops the run here.
    ch_samples = channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .map { row -> row instanceof List ? row[0] : row }

    STRUCTHOMOLOGY(ch_samples)
}
