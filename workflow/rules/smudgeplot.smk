#Run Smudgeplot
rule smudgeplot:
    output:
        "results/pre_assembly/smudgeplot/{tech}_smudgeplot.pdf"
    input:
        db_ktab = "results/pre_assembly/{tech}_reads.kmer.db.ktab",
        db_hist = "results/pre_assembly/{tech}_reads.kmer.db.hist"
    container:
        #"docker://quay.io/biocontainers/smudgeplot:0.5.3--py314h577a1d6_0"
        #this container below contains both smudgeplot and fastk since smudgeplot requires Logex which is distributed with FastK
        "docker://quay.io/sanger-tol/smudgeplot:0.5.4-c1"
    params:
        t = config["FastK_t"]
    resources:
        runtime = config["smudge_time"],
        mem_mb = config["smudge_mem"],
    shell:
        r"""
        smudgeplot hetmers results/pre_assembly/{wildcards.tech}_reads.kmer.db -L 10 -t {params.t} -o results/pre_assembly/smudgeplot/{wildcards.tech}_smudgeplot_pairs
        smudgeplot all results/pre_assembly/smudgeplot/{wildcards.tech}_smudgeplot_pairs_text.smu -o results/pre_assembly/smudgeplot/{wildcards.tech}
        """
