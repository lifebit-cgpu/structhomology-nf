/*
 * SOFTWARE_VERSIONS
 * Collate the versions.yml written by every process into one YAML file and
 * one TSV table, and add Nextflow and the pipeline version.
 */
process SOFTWARE_VERSIONS {
    label 'process_single'
    container 'quay.io/biocontainers/python:3.12.12'

    input:
    path(collated)
    val(nextflow_version)
    val(pipeline_version)
    path(script)

    output:
    path 'software_versions.yml', emit: yml
    path 'software_versions.tsv', emit: tsv

    script:
    """
    python3 ${script} \\
        --input ${collated} \\
        --nextflow ${nextflow_version} \\
        --pipeline ${pipeline_version} \\
        --process ${task.process} \\
        --out-prefix software_versions
    """

    stub:
    """
    touch software_versions.yml software_versions.tsv
    """
}
