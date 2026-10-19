
#Count k-mers for genomescope and smudgeplot
rule count_kmers:
    output:
        db_ktab = "results/pre_assembly/{tech}_reads.kmer.db.ktab",
        db_hist = "results/pre_assembly/{tech}_reads.kmer.db.hist",
        hist = "results/pre_assembly/{tech}_reads.histo",
    input:
        "input/{tech}.fastq.gz"
    threads: config["FastK_t"],
    resources:
        runtime = config["FastK_time"],
        mem_mb = config["FastK_m"],
    container:
        "docker://quay.io/biocontainers/fastk:1.2--h71df26d_1"
    params:
        k = config["FastK_k"],
        t = config["FastK_t"],
        m = config["FastK_m"]
    shell:
        r"""
        mkdir -p {wildcards.tech}_fastk_temp
        FastK -v -t1 \
         -T{params.t} \
         -k{params.k} \
         -M{params.m} \
         -P{wildcards.tech}_fastk_temp \
         {input} \
         -Nresults/pre_assembly/{wildcards.tech}_reads.kmer.db
        Histex -G results/pre_assembly/{wildcards.tech}_reads.kmer.db \
         > {output.hist}
        rm -r {wildcards.tech}_fastk_temp
        """
