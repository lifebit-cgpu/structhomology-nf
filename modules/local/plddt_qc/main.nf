/*
 * PLDDT_QC
 * Mean pLDDT, AlphaFold confidence-band fractions and the search gate for one
 * model. The gate decision (PASS / FAIL_MEAN / FAIL_FRAC / NO_MODEL) is emitted
 * as an environment value so the workflow can filter on it.
 */
process PLDDT_QC {
    tag "${meta.id}"
    label 'process_single'
    container 'quay.io/biocontainers/python:3.12.12'

    input:
    tuple val(meta), path(model), path(confidence), path(fasta), path(record)
    path(script)

    output:
    tuple val(meta), env('QC_FLAG')                    , emit: flag
    tuple val(meta), path("${meta.id}.plddt_qc.tsv")  , emit: tsv
    tuple val(meta), path("${meta.id}.plddt_qc.json") , emit: json
    path 'versions.yml'                               , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    python3 ${script} \\
        --accession ${meta.id} \\
        --group ${meta.group ?: 'query'} \\
        --model ${model} \\
        --confidence ${confidence} \\
        --record ${record} \\
        --out-prefix ${meta.id} \\
        ${args}

    QC_FLAG=\$(cat ${meta.id}.qc_flag.txt)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c 'import platform; print(platform.python_version())')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.plddt_qc.tsv
    echo '{}' > ${meta.id}.plddt_qc.json
    QC_FLAG=PASS
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
