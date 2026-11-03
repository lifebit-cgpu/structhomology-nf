/*
 * STRUCTHOMOLOGY: main analysis workflow
 *
 *  1. AFDB_FETCH          AlphaFold DB record + model; version resolved from the API
 *  2. PLDDT_QC            mean pLDDT, confidence-band fractions, search gate
 *  3. target database     one of: --foldseek_db (local copy), --target_pdb_ids
 *                         (small DB built from RCSB files), or FOLDSEEK_DATABASES
 *                         (downloads --foldseek_db_name after checking `foldseek databases`)
 *  4. DB_TO_FASTA         amino-acid sequences of exactly the same targets
 *  5. FOLDSEEK_EASYSEARCH structure arm;  MMSEQS_EASYSEARCH sequence arm
 *  6. ANNOTATE_HITS       merge arms, flag shared / structure-only / sequence-only, annotate
 *  7. SUMMARISE           per-protein summary + combined hit table
 *  8. CHECK_CONTROLS      positive-control acceptance test (optional, --controls)
 *  9. SOFTWARE_VERSIONS   versions.yml from every process -> one table
 */

include { AFDB_FETCH          } from '../modules/local/afdb_fetch/main'
include { PLDDT_QC            } from '../modules/local/plddt_qc/main'
include { FOLDSEEK_DATABASES  } from '../modules/local/foldseek_databases/main'
include { FETCH_PDB_TARGETS   } from '../modules/local/fetch_pdb_targets/main'
include { DB_TO_FASTA         } from '../modules/local/db_to_fasta/main'
include { ANNOTATE_HITS       } from '../modules/local/annotate_hits/main'
include { SUMMARISE           } from '../modules/local/summarise/main'
include { CHECK_CONTROLS      } from '../modules/local/check_controls/main'
include { SOFTWARE_VERSIONS   } from '../modules/local/software_versions/main'
include { FOLDSEEK_CREATEDB   } from '../modules/nf-core/foldseek/createdb/main'
include { FOLDSEEK_EASYSEARCH } from '../modules/nf-core/foldseek/easysearch/main'
include { MMSEQS_CREATEDB     } from '../modules/nf-core/mmseqs/createdb/main'
include { MMSEQS_EASYSEARCH   } from '../modules/nf-core/mmseqs/easysearch/main'

workflow STRUCTHOMOLOGY {

    take:
    ch_samples // channel: val(meta), meta = [ id: UniProt accession, group: 'query' | 'control' ]

    main:
    ch_versions = channel.empty()
    def bin = "${projectDir}/bin"

    // 1. AlphaFold DB record, model, confidence JSON and sequence
    AFDB_FETCH(ch_samples, file("${bin}/afdb_fetch.py", checkIfExists: true))
    ch_versions = ch_versions.mix(AFDB_FETCH.out.versions)

    // 2. pLDDT summary and gate
    PLDDT_QC(AFDB_FETCH.out.model, file("${bin}/plddt_qc.py", checkIfExists: true))
    ch_versions = ch_versions.mix(PLDDT_QC.out.versions)

    ch_pass = AFDB_FETCH.out.model
        .join(PLDDT_QC.out.flag)
        .filter { _meta, _model, _conf, _fasta, _record, flag -> flag == 'PASS' }
        .map { meta, model, _conf, fasta, _record, _flag -> [ meta, model, fasta ] }

    // 3. Target structure database
    def ch_db = channel.empty()
    if (params.foldseek_db) {
        def db_dir = file(params.foldseek_db, checkIfExists: true)
        ch_db = channel.value([ [ id: db_dir.name, name: db_dir.name, source: 'user-supplied directory' ], db_dir ])
    }
    else if (params.target_pdb_ids) {
        def ids_file = file(params.target_pdb_ids, checkIfExists: true)
        FETCH_PDB_TARGETS(
            [ id: 'custom_targets', name: ids_file.name, source: 'RCSB mmCIF download' ],
            ids_file,
            file("${bin}/fetch_pdb_targets.py", checkIfExists: true)
        )
        FOLDSEEK_CREATEDB(FETCH_PDB_TARGETS.out.structures)
        ch_db = FOLDSEEK_CREATEDB.out.db
        ch_versions = ch_versions.mix(FETCH_PDB_TARGETS.out.versions, FOLDSEEK_CREATEDB.out.versions)
    }
    else {
        def db_name = params.foldseek_db_name.toString()
        FOLDSEEK_DATABASES([ id: db_name.replaceAll(/[^A-Za-z0-9]+/, '_').toLowerCase(), name: db_name, source: 'foldseek databases' ])
        ch_db = FOLDSEEK_DATABASES.out.db
        FOLDSEEK_DATABASES.out.list.collectFile(name: 'foldseek_databases_list.txt', storeDir: "${params.outdir}/database")
        ch_versions = ch_versions.mix(FOLDSEEK_DATABASES.out.versions)
    }

    // 4. The same targets as amino-acid sequences, so both arms search one set
    DB_TO_FASTA(ch_db)
    MMSEQS_CREATEDB(DB_TO_FASTA.out.fasta)
    ch_seqdb = MMSEQS_CREATEDB.out.db
    ch_versions = ch_versions.mix(DB_TO_FASTA.out.versions, MMSEQS_CREATEDB.out.versions)

    // 5. Structure arm and sequence arm
    FOLDSEEK_EASYSEARCH(ch_pass.map { meta, model, _fasta -> [ meta, model ] }, ch_db)
    MMSEQS_EASYSEARCH(ch_pass.map { meta, _model, fasta -> [ meta, fasta ] }, ch_seqdb)
    ch_versions = ch_versions.mix(FOLDSEEK_EASYSEARCH.out.versions, MMSEQS_EASYSEARCH.out.versions)

    // 6. Merge arms, classify, annotate
    ch_merge = FOLDSEEK_EASYSEARCH.out.aln
        .join(MMSEQS_EASYSEARCH.out.tsv)
        .join(PLDDT_QC.out.json)
    ANNOTATE_HITS(ch_merge, file("${bin}/annotate_hits.py", checkIfExists: true))
    ch_versions = ch_versions.mix(ANNOTATE_HITS.out.versions)

    // 7. Run-level tables
    SUMMARISE(
        ANNOTATE_HITS.out.top_hits.map { _meta, f -> f }.collect().ifEmpty([]),
        ANNOTATE_HITS.out.summary.map { _meta, f -> f }.collect().ifEmpty([]),
        PLDDT_QC.out.tsv.map { _meta, f -> f }.collect(),
        DB_TO_FASTA.out.info,
        file("${bin}/summarise.py", checkIfExists: true)
    )
    ch_versions = ch_versions.mix(SUMMARISE.out.versions)

    // 8. Positive controls: the run fails if a control misses its known relatives
    if (params.controls) {
        CHECK_CONTROLS(
            SUMMARISE.out.combined,
            ANNOTATE_HITS.out.all_hits.map { _meta, f -> f }.collect().ifEmpty([]),
            file(params.controls, checkIfExists: true),
            file("${bin}/check_controls.py", checkIfExists: true)
        )
        ch_versions = ch_versions.mix(CHECK_CONTROLS.out.versions)
    }

    // 9. One software versions table
    SOFTWARE_VERSIONS(
        ch_versions.collectFile(name: 'collated_versions.yml', newLine: true, sort: true),
        workflow.nextflow.version.toString(),
        workflow.manifest.version ?: 'dev',
        file("${bin}/software_versions.py", checkIfExists: true)
    )

    emit:
    combined = SUMMARISE.out.combined
    summary  = SUMMARISE.out.per_protein
    versions = SOFTWARE_VERSIONS.out.yml
}
