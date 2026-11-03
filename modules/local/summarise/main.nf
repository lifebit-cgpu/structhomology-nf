/*
 * SUMMARISE
 * Run-level tables: one combined hit table, one per-protein summary (every
 * protein, including those stopped by the pLDDT gate), the pLDDT QC table and
 * a JSON run summary with the database details.
 */
process SUMMARISE {
    label 'process_single'
    container 'quay.io/biocontainers/python:3.12.12'

    input:
    path(top_hits , stageAs: 'top_hits/*')
    path(summaries, stageAs: 'summaries/*')
    path(qc       , stageAs: 'qc/*')
    path(db_info)
    path(script)

    output:
    path 'combined_hits.tsv'      , emit: combined
    path 'per_protein_summary.tsv', emit: per_protein
    path 'plddt_qc_summary.tsv'   , emit: qc
    path 'run_summary.json'       , emit: run
    path 'versions.yml'           , emit: versions

    script:
    """
    mkdir -p top_hits summaries qc
    python3 ${script} \\
        --top-hits-dir top_hits \\
        --summaries-dir summaries \\
        --qc-dir qc \\
        --db-info ${db_info} \\
        --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c 'import platform; print(platform.python_version())')
    END_VERSIONS
    """

    stub:
    """
    touch combined_hits.tsv per_protein_summary.tsv plddt_qc_summary.tsv
    echo '{}' > run_summary.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
