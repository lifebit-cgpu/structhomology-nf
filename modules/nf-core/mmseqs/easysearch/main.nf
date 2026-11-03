/*
 * MMSEQS_EASYSEARCH
 * Adapted from nf-core/modules (modules/nf-core/mmseqs/easysearch), same
 * container family and interface: tuple(meta, query fasta), tuple(meta_db,
 * db directory) -> tuple(meta, <id>.m8). The database prefix is detected from
 * its .lookup file.
 */
process MMSEQS_EASYSEARCH {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/mmseqs2:18.8cc5c--hd6d6fdc_0'

    input:
    tuple val(meta), path(query)
    tuple val(meta_db), path(db)

    output:
    tuple val(meta), path("${meta.id}.m8"), emit: tsv
    path 'versions.yml'                   , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    DB=\$(find -L ${db}/ -maxdepth 1 -name '*.lookup' | head -n 1 | sed 's/\\.lookup\$//')
    if [ -z "\$DB" ]; then
        echo "ERROR: no MMseqs2 database (*.lookup) found in ${db}" >&2
        exit 1
    fi

    mmseqs \\
        easy-search \\
        ${query} \\
        \$DB \\
        ${meta.id}.m8 \\
        tmpFolder \\
        --threads ${task.cpus} \\
        ${args}
    rm -rf tmpFolder

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs version)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.m8
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: stub
    END_VERSIONS
    """
}
