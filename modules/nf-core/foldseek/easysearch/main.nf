/*
 * FOLDSEEK_EASYSEARCH
 * Adapted from nf-core/modules (modules/nf-core/foldseek/easysearch), same
 * container and interface: tuple(meta, structure), tuple(meta_db, db directory)
 * -> tuple(meta, <id>.m8). Change: the database prefix is detected from the
 * *_ss.dbtype file so databases downloaded by `foldseek databases` (prefix
 * differs from meta_db.id) also work.
 */
process FOLDSEEK_EASYSEARCH {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/foldseek:10.941cd33--h5021889_1'

    input:
    tuple val(meta), path(pdb)
    tuple val(meta_db), path(db)

    output:
    tuple val(meta), path("${meta.id}.m8"), emit: aln
    path 'versions.yml'                   , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    DB=\$(find -L ${db}/ -maxdepth 1 -name '*_ss.dbtype' | head -n 1 | sed 's/_ss\\.dbtype\$//')
    if [ -z "\$DB" ]; then
        echo "ERROR: no Foldseek database (*_ss.dbtype) found in ${db}" >&2
        exit 1
    fi

    foldseek \\
        easy-search \\
        ${pdb} \\
        \$DB \\
        ${meta.id}.m8 \\
        tmpFolder \\
        --threads ${task.cpus} \\
        ${args}
    rm -rf tmpFolder

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        foldseek: \$(foldseek version)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.m8
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        foldseek: stub
    END_VERSIONS
    """
}
