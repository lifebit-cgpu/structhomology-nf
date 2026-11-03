/*
 * FOLDSEEK_EASYSEARCH
 * Adapted from nf-core/modules (modules/nf-core/foldseek/easysearch), same
 * container and interface: tuple(meta, structure), tuple(meta_db, db directory)
 * -> tuple(meta, <id>.m8). Change: the database prefix is detected inside the
 * directory so databases downloaded by `foldseek databases` (prefix differs
 * from meta_db.id) also work. If several candidate prefixes exist, the
 * shortest one that is a complete Foldseek structure database (<p>.dbtype,
 * <p>_ss.dbtype, <p>_ca.dbtype and <p>.lookup) is used, and the choice is
 * logged to stderr.
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
    DB=""
    for f in \$(ls ${db}/*_ss.dbtype 2>/dev/null | awk '{ print length(\$0), \$0 }' | sort -n | cut -d' ' -f2-); do
        p=\${f%_ss.dbtype}
        echo "candidate Foldseek prefix: \$p" >&2
        if [ -z "\$DB" ] && [ -e "\${p}.dbtype" ] && [ -e "\${p}_ca.dbtype" ] && [ -e "\${p}.lookup" ]; then DB=\$p; fi
    done
    if [ -z "\$DB" ]; then
        echo "ERROR: no complete Foldseek structure database (<p>.dbtype, <p>_ss.dbtype, <p>_ca.dbtype, <p>.lookup) found in ${db}" >&2
        exit 1
    fi
    echo "using Foldseek database prefix: \$DB" >&2

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
