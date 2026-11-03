/*
 * FOLDSEEK_DATABASES
 * Runs `foldseek databases` to list what is available, stops if the requested
 * name is not on that list, then downloads it. Cached with storeDir (see
 * conf/modules.config) so later runs reuse the copy.
 * nf-core/modules has no foldseek/databases module; this follows the same
 * conventions as mmseqs/databases.
 */
process FOLDSEEK_DATABASES {
    tag "${meta.name}"
    label 'process_medium'
    container 'quay.io/biocontainers/foldseek:10.941cd33--h5021889_1'

    input:
    val(meta) // [ id: safe directory name, name: name as listed by `foldseek databases` ]

    output:
    tuple val(meta), path("${meta.id}")              , emit: db
    path("${meta.id}.databases_list.txt")            , emit: list
    path("${meta.id}.versions.yml")                  , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    foldseek databases > ${meta.id}.databases_list.txt 2>&1 || true

    if ! grep -qwF -- "${meta.name}" ${meta.id}.databases_list.txt; then
        echo "ERROR: '${meta.name}' is not offered by 'foldseek databases'. Available databases:" >&2
        cat ${meta.id}.databases_list.txt >&2
        exit 1
    fi

    mkdir -p ${meta.id} tmp_download
    foldseek databases \\
        ${meta.name} \\
        ${meta.id}/${meta.id} \\
        tmp_download \\
        --threads ${task.cpus} \\
        ${args}
    rm -rf tmp_download

    cat <<-END_VERSIONS > ${meta.id}.versions.yml
    "${task.process}":
        foldseek: \$(foldseek version)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}
    touch ${meta.id}/${meta.id} ${meta.id}.databases_list.txt
    cat <<-END_VERSIONS > ${meta.id}.versions.yml
    "${task.process}":
        foldseek: stub
    END_VERSIONS
    """
}
