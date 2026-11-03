/*
 * ANNOTATE_HITS
 * Merge the Foldseek and MMseqs2 hits of one protein by target, call
 * significance in each arm, classify every target as shared, structure-only,
 * sequence-only or weak, and attach functional annotation from RCSB PDB and
 * UniProt. Flags hits that are themselves poorly annotated and hits to the
 * query protein itself.
 */
process ANNOTATE_HITS {
    tag "${meta.id}"
    label 'process_single'
    container 'quay.io/biocontainers/python:3.12.12'

    input:
    tuple val(meta), path(foldseek_hits, stageAs: 'foldseek/*'), path(mmseqs_hits, stageAs: 'mmseqs/*'), path(qc_json)
    path(script)

    output:
    tuple val(meta), path("${meta.id}.top_hits.tsv"), emit: top_hits
    tuple val(meta), path("${meta.id}.all_hits.tsv"), emit: all_hits
    tuple val(meta), path("${meta.id}.summary.json"), emit: summary
    path 'versions.yml'                              , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    python3 ${script} \\
        --accession ${meta.id} \\
        --group ${meta.group ?: 'query'} \\
        --foldseek ${foldseek_hits} \\
        --mmseqs ${mmseqs_hits} \\
        --qc ${qc_json} \\
        --out-prefix ${meta.id} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c 'import platform; print(platform.python_version())')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.top_hits.tsv ${meta.id}.all_hits.tsv
    echo '{}' > ${meta.id}.summary.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
