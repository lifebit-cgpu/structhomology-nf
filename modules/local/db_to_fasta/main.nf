/*
 * DB_TO_FASTA
 * Extract the amino-acid sequences of every entry in the Foldseek target
 * database (same identifiers), so the MMseqs2 arm searches exactly the same
 * targets. Also records the database version, size on disk and entry count.
 * Prefix detection matches FOLDSEEK_EASYSEARCH: the shortest complete
 * Foldseek structure database in the directory. Size on disk counts only the
 * files belonging to that database.
 */
process DB_TO_FASTA {
    tag "${meta.id}"
    label 'process_low'
    container 'quay.io/biocontainers/foldseek:10.941cd33--h5021889_1'

    input:
    tuple val(meta), path(db)

    output:
    tuple val(meta), path("${meta.id}.targets.fasta"), emit: fasta
    path("${meta.id}.db_info.json")                   , emit: info
    path 'versions.yml'                               , emit: versions

    script:
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

    foldseek convert2fasta \$DB ${meta.id}.targets.fasta

    SIZE=\$(for f in \$DB \$DB.* \${DB}_ss \${DB}_ss.* \${DB}_ca \${DB}_ca.* \${DB}_h \${DB}_h.* \${DB}_clu \${DB}_clu.* \${DB}_mapping \${DB}_taxonomy; do [ -f "\$f" ] && stat -L -c %s "\$f"; done | awk '{s+=\$1} END {printf "%d", s}')
    NENT=\$(wc -l < \$DB.lookup | tr -d ' ')
    NSEQ=\$(grep -c '^>' ${meta.id}.targets.fasta || true)
    VER=\$(for f in \$DB.version \${DB}_version; do [ -f "\$f" ] && cat "\$f"; done | tr '\\n\\t"' '   ' | sed 's/ *\$//')
    if [ -z "\$VER" ]; then
        VER=\$(find -L ${db}/ -maxdepth 1 -type f -iname '*version*' -exec cat {} + 2>/dev/null | tr '\\n\\t"' '   ' | sed 's/ *\$//')
    fi
    FSV=\$(foldseek version)

    cat <<-END_INFO > ${meta.id}.db_info.json
    {"name": "${meta.name ?: meta.id}", "id": "${meta.id}", "source": "${meta.source ?: ''}", "prefix": "\$(basename \$DB)", "version": "\$VER", "size_bytes": \$SIZE, "n_entries": \$NENT, "n_sequences": \$NSEQ, "foldseek_version": "\$FSV"}
    END_INFO

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        foldseek: \$FSV
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.targets.fasta
    echo '{}' > ${meta.id}.db_info.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        foldseek: stub
    END_VERSIONS
    """
}
