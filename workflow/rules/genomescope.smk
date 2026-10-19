#Run Genomescope
rule genomescope:
    output:
        "results/pre_assembly/genomescope/{tech}_genomescope_ploidy{ploidy}.out"
    input:
        #rules.count_kmers.output.hist
        histo = "results/pre_assembly/{tech}_reads.histo",
    resources:
        runtime = config["GS_time"],
        mem_mb = config["GS_mem"],
    params:
        ploidy_range = config["GS_ploidy_range"],
        k = config["FastK_k"]
    container:
        "docker://quay.io/biocontainers/genomescope2:2.1.0--py313r44hdfd78af_0"
    shell:
        r"""
        genomescope2 \
         -i {input.histo} \
         -o results/pre_assembly/genomescope/{wildcards.tech}_ploidy{wildcards.ploidy} \
         -k {params.k} \
         -p {wildcards.ploidy} \
         1> results/pre_assembly/genomescope/{wildcards.tech}_genomescope_ploidy{wildcards.ploidy}.out \
         2> results/pre_assembly/genomescope/{wildcards.tech}_genomescope_ploidy{wildcards.ploidy}.err
        """
