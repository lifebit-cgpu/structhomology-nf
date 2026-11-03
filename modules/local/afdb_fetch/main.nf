/*
 * AFDB_FETCH
 * Fetch the AlphaFold DB record for one UniProt accession, resolve the model
 * version from the API ('latestVersion' of the canonical AF-<acc>-F1 record),
 * and download the model, the per-residue confidence JSON and the sequence.
 */
process AFDB_FETCH {
    tag "${meta.id}"
    label 'process_single'
    container 'quay.io/biocontainers/python:3.12.12'

    input:
    val(meta)
    path(script)

    output:
    tuple val(meta), path("${meta.id}.model.*"), path("${meta.id}.confidence.json"), path("${meta.id}.fasta"), path("${meta.id}.record.json"), emit: model
    path 'versions.yml', emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    python3 ${script} \\
        --accession ${meta.id} \\
        --out-prefix ${meta.id} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c 'import platform; print(platform.python_version())')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.model.pdb ${meta.id}.confidence.json ${meta.id}.fasta
    echo '{}' > ${meta.id}.record.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
