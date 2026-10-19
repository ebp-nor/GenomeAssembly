rule read_statistics:
    output:
        stats = "results/stats/{tech}.stats"
    input:
        reads = "input/{tech}.fastq.gz"
    container:
        "docker://quay.io/biocontainers/seqkit:2.14.0--hb192632_0"
    shell:
        r"""
        seqkit stats -a {input.reads} > {output.stats}
        """
