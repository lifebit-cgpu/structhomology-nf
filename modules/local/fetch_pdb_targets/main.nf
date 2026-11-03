/*
 * FETCH_PDB_TARGETS
 * Download mmCIF files for a list of PDB IDs from RCSB, to build a small
 * Foldseek target database (used by -profile test).
 */
process FETCH_PDB_TARGETS {
    tag "${meta.id}"
    label 'process_single'
    container 'quay.io/biocontainers/python:3.12.12'

    input:
    val(meta)
    path(ids)
    path(script)

    output:
    tuple val(meta), path('targets'), emit: structures
    path 'targets_manifest.tsv'     , emit: manifest
    path 'versions.yml'             , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    python3 ${script} \\
        --ids ${ids} \\
        --outdir targets \\
        --manifest targets_manifest.tsv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c 'import platform; print(platform.python_version())')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p targets
    touch targets/1abc.cif.gz targets_manifest.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
