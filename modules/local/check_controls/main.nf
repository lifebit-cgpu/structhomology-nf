/*
 * CHECK_CONTROLS
 * Positive-control acceptance test. A control passes when one of its top 10
 * Foldseek hits is a known relative (expected PDB ID or family name) AND is
 * significant (E <= --sig_evalue or TM-score >= --min_tmscore). With
 * --fail_on_control_failure the run stops if any control fails.
 */
process CHECK_CONTROLS {
    label 'process_single'
    container 'quay.io/biocontainers/python:3.12.12'

    input:
    path(combined)
    path(all_hits, stageAs: 'all_hits/*')
    path(controls)
    path(script)

    output:
    path 'control_check.tsv', emit: tsv
    path 'versions.yml'     , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p all_hits
    python3 ${script} \\
        --combined ${combined} \\
        --all-hits-dir all_hits \\
        --controls ${controls} \\
        --out control_check.tsv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c 'import platform; print(platform.python_version())')
    END_VERSIONS
    """

    stub:
    """
    touch control_check.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
