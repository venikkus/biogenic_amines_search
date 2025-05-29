configfile: "config.yaml"

SOURCE = config["source"]
INPUT_TYPE = config["input_type"]
IDS = config["id"]  # Список ID
GENES = config['genes_of_interest']

rule all:
    input:
        expand("assembly/{sample}/contigs.fasta", sample=IDS),
        expand("quast_output/{sample}/report.tsv", sample=IDS),
        expand("prokka_output/{sample}/{sample}.gff", sample=IDS),
        expand("prokka_output/{sample}/{sample}.tsv", sample=IDS),
        expand("decarboxylases/amines_{sample}/amines_{sample}.tsv", sample=IDS),
        "summary/presence_absence_matrix.tsv",
        "summary/presence_absence_heatmap.png"

if INPUT_TYPE == "sra":
    if SOURCE == "offline":
        rule copy_offline_reads:
            input:
                r1 = config["offline_reads"][0],
                r2 = config["offline_reads"][1]
            output:
                "input_reads/{sample}_1.fastq",
                "input_reads/{sample}_2.fastq"
            shell:
                """
                mkdir -p input_reads
                cp {input.r1} input_reads/{wildcards.sample}_1.fastq
                cp {input.r2} input_reads/{wildcards.sample}_2.fastq
                """

    if SOURCE == "online":
        rule download_sra:
            output:
                "sra/{sample}/{sample}.sra"
            shell:
                "prefetch {wildcards.sample} -O sra/"

        rule split_files:
            input:
                "sra/{sample}/{sample}.sra"
            output:
                "input_reads/{sample}_1.fastq",
                "input_reads/{sample}_2.fastq"
            shell:
                "fasterq-dump --split-files {input} --outdir input_reads/"

    rule trimming_reads:
        input:
            "input_reads/{sample}_1.fastq",
            "input_reads/{sample}_2.fastq"
        output:
            "clean_data/{sample}_1P.fastq",
            "clean_data/{sample}_2P.fastq"
        shell:
            """
            mkdir -p clean_data
            trimmomatic PE {input[0]} {input[1]} \
            clean_data/{wildcards.sample}_1P.fastq clean_data/{wildcards.sample}_1U.fastq \
            clean_data/{wildcards.sample}_2P.fastq clean_data/{wildcards.sample}_2U.fastq \
            LEADING:20 TRAILING:20 MINLEN:20
            """

    rule genome_assembly:
        input:
            "clean_data/{sample}_1P.fastq",
            "clean_data/{sample}_2P.fastq"
        output:
            "assembly/{sample}/contigs.fasta"
        shell:
            """
            mkdir -p assembly/{wildcards.sample}
            spades.py -1 {input[0]} -2 {input[1]} -o assembly/{wildcards.sample}
            """

    rule run_quast:
        input:
            "assembly/{sample}/contigs.fasta"
        output:
            directory("quast_output/{sample}")
        shell:
            """
            mkdir -p quast_output/{wildcards.sample}
            quast {input} -o quast_output/{wildcards.sample}
            """

if INPUT_TYPE == "assembly":
    if SOURCE == "offline":
        rule copy_offline_assembly:
            input:
                config["offline_assembly"]
            output:
                "assembly/{sample}/contigs.fasta"
            shell:
                """
                mkdir -p assembly/{wildcards.sample}
                cp {input} assembly/{wildcards.sample}/contigs.fasta
                """

    if SOURCE == "online":
        rule download_assembly:
            output:
                "assembly/{sample}/contigs.fasta"
            shell:
                """
                mkdir -p assembly/{wildcards.sample}
                datasets download genome accession {wildcards.sample} --filename {wildcards.sample}.zip
                unzip -o {wildcards.sample}.zip -d assembly/{wildcards.sample}
                # Ищем fasta файл (например, *.fna) в папке разархивированных данных
                find assembly/{wildcards.sample} -name '*.fna' -exec cp {{}} assembly/{wildcards.sample}/contigs.fasta \\;
                rm {wildcards.sample}.zip
                """

# общие
rule run_prokka:
    input:
        "assembly/{sample}/contigs.fasta"
    output:
        gff = "prokka_output/{sample}/{sample}.gff",
        tsv = "prokka_output/{sample}/{sample}.tsv",
    shell:
        """
        mkdir -p prokka_output/{wildcards.sample}
        prokka --force --outdir prokka_output/{wildcards.sample} --prefix {wildcards.sample} {input}
        """

rule amine_search:
    input:
        "prokka_output/{sample}/{sample}.tsv"
    output:
        "decarboxylases/amines_{sample}/amines_{sample}.tsv"
    shell:
        r"""
        mkdir -p decarboxylases/amines_{wildcards.sample}
        awk -F '\t' '$4 ~ /tdc|hdc|ldc|cad|odc|adi|puu|pot|trp|spe/ {{print $1 "\t" $4 "\t" $7}}' {input} > {output}
        """

rule presence_absence_matrix:
    input:
        sample_files=expand("decarboxylases/amines_{sample}/amines_{sample}.tsv", sample=IDS)
    output:
        "summary/presence_absence_matrix.tsv"
    run:
        import os
        import pandas as pd
        import re
        from collections import defaultdict
        
        os.makedirs("summary", exist_ok=True)
        gene_pattern = re.compile(r'tdc|hdc|ldc|cad|odc|adi|puu|pot|trp|spe')
        base_gene_names = set()

        for file in input.sample_files:
            if os.path.exists(file) and os.path.getsize(file) > 0:
                try:
                    df = pd.read_csv(file, sep='\t', header=None, names=["locus_tag", "gene", "product"])
                    # extracting base gene names (without suffixes)
                    base_names = [re.sub(r'_\d+$', '', gene) for gene in df["gene"] if gene_pattern.search(gene)]
                    base_gene_names.update(base_names)
                except (pd.errors.EmptyDataError, FileNotFoundError):
                    continue
        
        base_gene_names = sorted(base_gene_names)
        data = defaultdict(lambda: defaultdict(int))
        
        # copies count
        for file in input.sample_files:
            sample = os.path.basename(os.path.dirname(file)).replace("amines_", "")
            
            for gene in base_gene_names:
                data[sample][gene] = 0
            
            if os.path.exists(file) and os.path.getsize(file) > 0:
                try:
                    df = pd.read_csv(file, sep='\t', header=None, names=["locus_tag", "gene", "product"])
                    filtered_df = df[df["gene"].str.contains(gene_pattern, na=False)]
                    
                    for _, row in filtered_df.iterrows():
                        base_name = re.sub(r'_\d+$', '', row["gene"])
                        if base_name in base_gene_names:
                            data[sample][base_name] += 1
                except (pd.errors.EmptyDataError, FileNotFoundError):
                    continue
        
        result_df = pd.DataFrame.from_dict({k: dict(v) for k, v in data.items()}, orient='index')
        result_df = result_df.reindex(columns=base_gene_names)
        result_df.to_csv(output[0], sep='\t', index=True)

rule presence_absence_heatmap:
    input:
        tsv="summary/presence_absence_matrix.tsv",
    output:
        "summary/presence_absence_heatmap.png"
    shell:
        "python scripts/plot_heatmap.py {input.tsv} {output}"