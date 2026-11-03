/*
 * FOLDSEEK_CREATEDB
 * Adapted from nf-core/modules (modules/nf-core/foldseek/createdb), same
 * container and interface: tuple(meta, structures) -> tuple(meta, db directory),
 * the database prefix inside the directory being meta.id.
 */
process FOLDSEEK_CREATEDB {
    tag "${meta.id}"
    label 'process_low'
    container 'quay.io/biocontainers/foldseek:10.941cd33--h5021889_1'

    input:
    tuple val(meta), path(pdb)

    output:
    tuple val(meta), path("${meta.id}"), emit: db
    path 'versions.yml'                , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p ${meta.id}
    foldseek \\
        createdb \\
        ${pdb} \\
        ${meta.id}/${meta.id} \\
        --threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        foldseek: \$(foldseek version)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}
    touch ${meta.id}/${meta.id} ${meta.id}/${meta.id}_ss.dbtype
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        foldseek: stub
    END_VERSIONS
    """
}
