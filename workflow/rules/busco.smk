#Run busco on assembly
rule busco:
    output:
        directory(expand("results/{{dir}}/busco_{{assembly}}_{lineage}", lineage=config["busco_lineage"]))
    input:
        assembly = "results/{dir}/{assembly}.fa",
        busco = expand("{busco_db_dir}/{lineage}_odb{version}", busco_db_dir=config["busco_db_dir"], lineage=config["busco_lineage"], version=config["busco_db_version"])
    threads: config["busco_threads"],
    resources:
        mem_mb = config["busco_mem"],
    params:
        lineage = config["busco_lineage"],
        threads = config["busco_threads"]
    container:
        "docker://quay.io/biocontainers/busco:6.1.0--pyhdfd78af_2"
    shell:
        r"""
        busco \
            -i {input.assembly} \
            -l {input.busco} \
            -c {params.threads} \
            -m genome \
            --offline \
            --tar \
            --out_path results/{wildcards.dir} \
            -o busco_{wildcards.assembly}_{params.lineage} \
            --download_path /tmp
        """
