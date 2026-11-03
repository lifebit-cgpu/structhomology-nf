/*
 * DB_TO_FASTA
 * Extract the amino-acid sequences of every entry in the Foldseek target
 * database (same identifiers), so the MMseqs2 arm searches exactly the same
 * targets. Also records the database version, size on disk and entry count.
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
    DB=\$(find -L ${db}/ -maxdepth 1 -name '*_ss.dbtype' | head -n 1 | sed 's/_ss\\.dbtype\$//')
    if [ -z "\$DB" ]; then
        echo "ERROR: no Foldseek database (*_ss.dbtype) found in ${db}" >&2
        exit 1
    fi

    foldseek convert2fasta \$DB ${meta.id}.targets.fasta

    SIZE=\$(find -L ${db}/ -type f -exec stat -L -c %s {} + | awk '{s+=\$1} END {printf "%d", s}')
    NENT=\$(wc -l < \$DB.lookup | tr -d ' ')
    NSEQ=\$(grep -c '^>' ${meta.id}.targets.fasta || true)
    VER=\$(find -L ${db}/ -maxdepth 1 -type f -iname '*version*' -exec cat {} + 2>/dev/null | tr '\\n\\t"' '   ' | sed 's/ *\$//')
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
