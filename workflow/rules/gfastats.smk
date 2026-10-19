rule gfastats:
    output:
        assembly_stats = "results/{dir}/{assembly}.fa.stats",
    input:
        assembly = "results/{dir}/{assembly}.fa",
    container:
        "docker://quay.io/biocontainers/gfastats:1.3.11--h077b44d_0"
    shell:
        r"""
        gfastats {input.assembly} > {output.assembly_stats}
        """
