/*
 * MMSEQS_CREATEDB
 * Adapted from nf-core/modules (modules/nf-core/mmseqs/createdb), same
 * container family and interface: tuple(meta, fasta) -> tuple(meta, db directory).
 */
process MMSEQS_CREATEDB {
    tag "${meta.id}"
    label 'process_low'
    container 'quay.io/biocontainers/mmseqs2:18.8cc5c--hd6d6fdc_0'

    input:
    tuple val(meta), path(sequence)

    output:
    tuple val(meta), path("${meta.id}_seqdb"), emit: db
    path 'versions.yml'                      , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p ${meta.id}_seqdb
    mmseqs \\
        createdb \\
        ${sequence} \\
        ${meta.id}_seqdb/${meta.id} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs version)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}_seqdb
    touch ${meta.id}_seqdb/${meta.id} ${meta.id}_seqdb/${meta.id}.lookup
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: stub
    END_VERSIONS
    """
}
