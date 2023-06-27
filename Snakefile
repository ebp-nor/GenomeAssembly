#Link to the config.yaml containing all params and other information
configfile: "EBPNor_GenomeAssembly_profile/EBPNor_GenomeAssembly_config.yaml"

#Get working directory from config
species_dir = config["species_dir"]
if not species_dir.endswith("/"):
    species_dir += "/"

#Fetch and store all *.fastq.gz files in the input directory
FASTQ, = glob_wildcards(species_dir+"genomic_data/pacbio/{sample,[^/]+}.fastq.gz")
HIC, = glob_wildcards(species_dir+"genomic_data/hic/{sample}_1.fastq.gz")
HAPS = ["hap1", "hap2"]

#Define rules that do not need to be submitted to the cluster
localrules: concatenate_filtered_reads

#Summary rule to run full pipeline (pre-assembly, assembly, scaffolding)
rule all:
    input:
        expand(species_dir+"results/pre_assembly/genomescope/genomescope_ploidy{ploidy}.out", ploidy=config["GS_ploidy_range"]),
        expand(species_dir+"results/pre_assembly/smudgeplot/ploidy{ploidy}_smudgeplot.png", ploidy=config["smudge_ploidy"]),
        expand(species_dir+"results/assembly/busco_assembly.hic.{hap}.p_ctg_{lineage}", lineage=config["busco_lineage"], hap=HAPS),
        expand(species_dir+"results/scaffolding/busco_{hap}_scaffolds_final_{lineage}", lineage=config["busco_lineage"], hap=HAPS)

#Summary rule to run all pre-assembly steps
rule pre_assembly:
    input:
        expand(species_dir+"results/pre_assembly/genomescope/genomescope_ploidy{ploidy}.out", ploidy=config["GS_ploidy_range"]),
        expand(species_dir+"results/pre_assembly/smudgeplot/ploidy{ploidy}_smudgeplot.png", ploidy=config["smudge_ploidy"])

#General rule for running only primary assembly and busco
rule assembly:
    input:
        expand(species_dir+"results/assembly/busco_assembly.hic.{hap}.p_ctg_{lineage}", lineage=config["busco_lineage"], hap=HAPS)
        
#General rule for generated the scaffolded haplotype assemblies (withouth pre-assembly) and busco
rule scaffolding:
    input:
        expand(species_dir+"results/assembly/busco_assembly.hic.{hap}.p_ctg_{lineage}", lineage=config["busco_lineage"], hap=HAPS),
        expand(species_dir+"results/scaffolding/busco_{hap}_scaffolds_final_{lineage}", lineage=config["busco_lineage"], hap=HAPS)

#Filter PacBio HiFi reads
rule hifiadapterfilt:
    output:
        expand(species_dir+"genomic_data/pacbio/filtered/{sample}.filt.fastq.gz", sample=FASTQ)
    input:
        expand(species_dir+"genomic_data/pacbio/{sample}.fastq.gz", sample=FASTQ)
    resources:
        ntasks = config["adapterfilt_threads"]
#    conda:
#        config["adapterfilt_conda_env"]
    params:
        threads = config["adapterfilt_threads"],
        installdir = config["adapterfilt_install_dir"],
        modules = config["adapterfilt_modules"],
        workdir = species_dir,
        tempdir = species_dir+"temp"
    shell:
        r"""
        module load {params.modules}
        export PATH={params.installdir}/DB:{params.installdir}:$PATH
        mkdir -p {params.tempdir}
        export TMPDIR={params.tempdir}
        cd {params.workdir}/genomic_data/pacbio/
        pbadapterfilt.sh -t {params.threads} -o filtered
        rm -r {params.tempdir}
        """

#Count k-mers for genomescope and smudgeplot
rule count_kmers:
    output:
        species_dir+"results/pre_assembly/reads.histo"
    input:
        expand(species_dir+"genomic_data/pacbio/{sample}.fastq.gz", sample=FASTQ)
 #   conda:
 #       config["GS_conda_env"]
    resources:
        ntasks = config["KMC_t"]
    params:
        outdir = species_dir+"results/pre_assembly",
        k = config["KMC_k"],
        t = config["KMC_t"],
        m = config["KMC_m"],
        ci = config["KMC_ci"],
        cx = config["KMC_cx"],
        cs = config["KMC_cs"]
    shell:
        r"""
        mkdir -p {params.outdir}/tmp
        ls {input} > FILES
        kmc -k{params.k} -t{params.t} -m{params.m} -ci{params.ci} -cs{params.cs} @FILES {params.outdir}/reads {params.outdir}/tmp
        kmc_tools transform {params.outdir}/reads histogram {params.outdir}/reads.histo -cx{params.cx}
        rm -r {params.outdir}/tmp FILES
        """

#Run Genomescope
rule genomescope:
    output:
        expand(species_dir+"results/pre_assembly/genomescope/genomescope_ploidy{ploidy}.out", ploidy=config["GS_ploidy_range"])
    input:
        species_dir+"results/pre_assembly/reads.histo"
  #  conda:
  #      config["GS_conda_env"]
    params:
        outdir = species_dir+"results/pre_assembly/genomescope",
        ploidy_range = config["GS_ploidy_range"],
        k = config["KMC_k"]
    shell:
        r"""
        for ploidy in {params.ploidy_range} 
        do genomescope2 -i {input} -o {params.outdir}/output_ploidy$ploidy -k {params.k} -p $ploidy 1> {params.outdir}/genomescope_ploidy${{ploidy}}.out 2> {params.outdir}/genomescope_ploidy${{ploidy}}.err
        done
        """

#Run Smudgeplot
rule smudgeplot:
    output:
        expand(species_dir+"results/pre_assembly/smudgeplot/ploidy{ploidy}_smudgeplot.png", ploidy=config["smudge_ploidy"]),
    input:
        species_dir+"results/pre_assembly/reads.histo"
  #  conda:
  #      config["smudge_conda_env"]
    params:
        indir = species_dir+"results/pre_assembly",
        outdir = species_dir+"results/pre_assembly/smudgeplot",
        ploidy = config["smudge_ploidy"],
        kmc_path = config["smudge_path"],
        module = config["smudge_module"],
        k = config["KMC_k"]
    shell:
        r"""
        L=$(smudgeplot.py cutoff {input} L)
        U=$(smudgeplot.py cutoff {input} U)
        export PATH={params.kmc_path}/bin:$PATH
        kmc_tools transform {params.indir}/reads -ci$L -cx$U reduce {params.outdir}/kmcdb_L"$L"_U"$U"
        smudge_pairs {params.outdir}/kmcdb_L"$L"_U"$U" {params.outdir}/kmcdb_L"$L"_U"$U"_coverages.tsv {params.outdir}/kmcdb_L"$L"_U"$U"_pairs.tsv > {params.outdir}/kmcdb_L"$L"_U"$U"_familysizes.tsv
        module unload {params.module}
        smudgeplot.py plot {params.outdir}/kmcdb_L"$L"_U"$U"_coverages.tsv -o {params.outdir}/ploidy{params.ploidy} -k {params.k}
        """

#Concatenate all filtered reads for assembly
rule concatenate_filtered_reads:
    output:
        species_dir+"genomic_data/pacbio/filtered/concat.filt.fastq.gz"
    input:
        expand(species_dir+"genomic_data/pacbio/filtered/{sample}.filt.fastq.gz", sample=FASTQ)
    shell:
        "cat {input} > {output}"

#Perform genome assembly with hifiasm using pacbio reads and hic reads
rule assembly_hifiams:
    output:
        hap1 = species_dir+"results/assembly/assembly.hic.hap1.p_ctg.fa",
        hap2 = species_dir+"results/assembly/assembly.hic.hap2.p_ctg.fa"
    input:
        pacbio = species_dir+"genomic_data/pacbio/filtered/concat.filt.fastq.gz",
        hicF = expand(species_dir+"genomic_data/hic/{hic}_1.fastq.gz", hic=HIC),
        hicR = expand(species_dir+"genomic_data/hic/{hic}_2.fastq.gz", hic=HIC)
    resources:
        time = config["hifiasm_cluster_time"],
        mem_per_cpu = config["hifiasm_mem_per_cpu"],
        ntasks = config["hifiasm_threads"]
  #  conda:
   #     config["hifiasm_conda_env"]
    params:
        outdir = species_dir+"results/assembly",
        threads = config["hifiasm_threads"],
    shell:
        r"""
        hifiasm -o {params.outdir}/assembly -t{params.threads} --h1 {input.hicF} --h2 {input.hicR} {input.pacbio} \
        1> {params.outdir}/hifiasm_"`date +\%y\%m\%d_\%H\%M\%S`".out 2> {params.outdir}/hifiasm_"`date +\%y\%m\%d_\%H\%M\%S`".err
        awk '/^S/{{print ">"$2"\n"$3}}' {params.outdir}/assembly.hic.hap1.p_ctg.gfa | fold > {output.hap1}
        awk '/^S/{{print ">"$2"\n"$3}}' {params.outdir}/assembly.hic.hap2.p_ctg.gfa | fold > {output.hap2}
        """

#Download busco_db before submitting busco job
#rule download_busco_db:
#    output:
#        temporary(directory(expand("tmp/busco/{lineage}_odb10", lineage=config["busco_lineage"])))
    #conda:
     #   config["busco_conda_env"]
#    params:
#        lineage = config["busco_lineage"]
#    shell:
#        r"""
#        busco --download {params.lineage}_odb10
#        mv busco_downloads/lineages/{params.lineage}_odb10 {output}
#        rm -r busco*
#        """

#Run busco on assembly
rule busco:
    output:
        directory(expand(species_dir+"results/{{dir}}/busco_{{assembly}}_{lineage}", lineage=config["busco_lineage"]))
    input:
        assembly = species_dir+"results/{dir}/{assembly}.fa",
        busco = expand("{db_dir}/{lineage}_odb10", db_dir=config["busco_db_dir"], lineage=config["busco_lineage"])
    resources:
        mem_per_cpu = config["busco_mem"],
        ntasks = config["busco_threads"]
    #conda:
    #    config["busco_conda_env"]
    params:
        lineage = config["busco_lineage"],
        threads = config["busco_threads"],
        speciesdir = species_dir
    shell:
        r"""
        busco -i {input.assembly} -l {input.busco} -c {params.threads} -m genome --offline --out_path {params.speciesdir}/results/{wildcards.dir} \
        -o busco_{wildcards.assembly}_{params.lineage} --download_path tmp/{input.assembly}_busco
        rm -r tmp/{input.assembly}_busco
        """

#Create Meryl Database
rule create_meryl_db:
    output:
        directory(species_dir+"results/scaffolding/meryl/assembly.hic.{hap}.p_ctg.meryl"),
    input:
        species_dir+"results/assembly/assembly.hic.{hap}.p_ctg.fa",
    #conda:
     #   config["meryl_conda_env"]
    resources:
        ntasks = config["meryl_threads"]
    params:
        k = config["meryl_k"],
        t = config["meryl_threads"],
        m = config["meryl_memory"]
    shell:
        "meryl k={params.k} threads={params.t} memory={params.m} count output {output} {input}"

#Run meryl for filtering Hi-C reads before scaffolding
rule meryl_filter:
    output:
        hicFwoH1 = species_dir+"results/scaffolding/meryl/hicF_hap2.fastq.gz",
        hicRwoH1 = species_dir+"results/scaffolding/meryl/hicR_hap2.fastq.gz",
        hicFwoH2 = species_dir+"results/scaffolding/meryl/hicF_hap1.fastq.gz",
        hicRwoH2 = species_dir+"results/scaffolding/meryl/hicR_hap1.fastq.gz"
    input:
        hicF = expand(species_dir+"genomic_data/hic/{hic}_1.fastq.gz", hic=HIC),
        hicR = expand(species_dir+"genomic_data/hic/{hic}_2.fastq.gz", hic=HIC),
        meryl_hap1 = species_dir+"results/scaffolding/meryl/assembly.hic.hap1.p_ctg.meryl",
        meryl_hap2 = species_dir+"results/scaffolding/meryl/assembly.hic.hap2.p_ctg.meryl"
    #conda:
    #    config["meryl_conda_env"]
    resources:
        ntasks = config["meryl_threads"]
    params:
        outdir = species_dir+"results/scaffolding/meryl"
    shell:
        r"""
        meryl difference {input.meryl_hap1} {input.meryl_hap2} output {params.outdir}/only_hap1.meryl
        meryl difference {input.meryl_hap2} {input.meryl_hap1} output {params.outdir}/only_hap2.meryl
        meryl-lookup -exclude -sequence {input.hicF} {input.hicR} -mers {params.outdir}/only_hap1.meryl -output {output.hicFwoH1} {output.hicRwoH1}
        meryl-lookup -exclude -sequence {input.hicF} {input.hicR} -mers {params.outdir}/only_hap2.meryl -output {output.hicFwoH2} {output.hicRwoH2}
        """

#Prepare filtered Hi-C reads for scaffolding
rule prepare_hic:
    output:
        species_dir+"results/scaffolding/{hap}_hic_markdup.sort_n.bam"
    input:
        assembly = species_dir+"results/assembly/assembly.hic.{hap}.p_ctg.fa",
        hicF = species_dir+"results/scaffolding/meryl/hicF_{hap}.fastq.gz",
        hicR = species_dir+"results/scaffolding/meryl/hicR_{hap}.fastq.gz"
    #conda:
    #    config["yahs_conda_env"]
    resources:
        ntasks = config["yahs_threads"],
    params:
        t = int(config["yahs_threads"]) - 6,
        outdir =species_dir+ "results/scaffolding"
    shell:
        r"""
        bwa index {input.assembly}
        samtools faidx {input.assembly} 
        bwa mem -t {params.t} -R '@RG\tSM:{wildcards.hap}\tID:{wildcards.hap}' -5SPM {input.assembly} {input.hicF} {input.hicR} \
        | samtools view -buS - \
        | samtools sort -@3 -n -T {params.outdir}/{wildcards.hap}_tmp_n -O bam - \
        | samtools fixmate -mr - - \
        | samtools sort -@3 -T {params.outdir}/{wildcards.hap}_hic_tmp -O bam \
        | samtools markdup -rsS - - 2> {params.outdir}/{wildcards.hap}_hic_markdup.stats \
        | samtools sort -n -@3 -n -T {params.outdir}/{wildcards.hap}_temp_n -O bam > {output}
        """

#Run scaffolding with yahs
rule scaffold_hap:
    output:
        species_dir+"results/scaffolding/{hap}_scaffolds_final.fa"
    input:
        assembly = species_dir+"results/assembly/assembly.hic.{hap}.p_ctg.fa",
        bamfile = species_dir+"results/scaffolding/{hap}_hic_markdup.sort_n.bam"
    #conda:
    #    config["yahs_conda_env"]
    resources:
        ntasks = config["yahs_threads"],
        time = config["yahs_cluster_time"]
    params:
        outdir = species_dir+"results/scaffolding"
    shell:
        r"""
        yahs {input.assembly} {input.bamfile} -o {params.outdir}/{wildcards.hap} 1> {params.outdir}/yahs_{wildcards.hap}.out 2> {params.outdir}/yahs_{wildcards.hap}.err
        gfastats {output} > {params.outdir}/{wildcards.hap}_scaffolds_final.stats
        """
