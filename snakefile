import os
configfile: "config.yaml"

ASSEMBLY_IDS = config["assembly_samples"]["id"]
SRA_IDS = config["sra_samples"]["id"]
ALL_IDS = ASSEMBLY_IDS + SRA_IDS

ASSEMBLY_SOURCE = config["assembly_samples"]["source"]
SRA_SOURCE = config["sra_samples"]["source"]
GENES_FILE = config["genes_of_interest"]

ASSEMBLY_LOCAL_DIR = config["assembly_samples"].get("local_path", "local_assemblies")
SRA_LOCAL_DIR = config["sra_samples"].get("local_path", "local_sra")

GENES = []
if os.path.exists(GENES_FILE):
    with open(GENES_FILE) as f:
        GENES = [line.strip() for line in f]

rule all:
    input:
        expand("assembly/{sample}/contigs.fasta", sample=ALL_IDS),
        expand("quast_output/{sample}/report.tsv", sample=ALL_IDS),
        expand("prokka_output/{sample}/{sample}.gff", sample=ALL_IDS),
        expand("prokka_output/{sample}/{sample}.tsv", sample=ALL_IDS),
        expand("decarboxylases/amines_{sample}/amines_{sample}.tsv", sample=ALL_IDS),
        "summary/presence_absence_matrix.tsv",
        "summary/presence_absence_heatmap.png",
        "reports/full_report.html"

if SRA_IDS:
    if SRA_SOURCE == "online":
        rule download_sra:
            output: 
                "sra/{sample}.sra"
            wildcard_constraints:
                sample="|".join(SRA_IDS)
            shell:
                "mkdir -p sra && "
                "prefetch {wildcards.sample} --output-file {output}"

        rule split_files:
            input:
                "sra/{sample}.sra"
            output:
                r1 = "input_reads/{sample}_1.fastq",
                r2 = "input_reads/{sample}_2.fastq"
            wildcard_constraints:
                sample="|".join(SRA_IDS)
            shell:
                "mkdir -p input_reads && "
                "fasterq-dump {input} -O input_reads/ --split-files"
    
    elif SRA_SOURCE == "offline":
        rule local_sra_input:
            input:
                r1 = os.path.join(SRA_LOCAL_DIR, "{sample}/{sample}_1.fastq"),
                r2 = os.path.join(SRA_LOCAL_DIR, "{sample}/{sample}_2.fastq")
            output:
                r1 = "input_reads/{sample}_1.fastq",
                r2 = "input_reads/{sample}_2.fastq"
            wildcard_constraints:
                sample="|".join(SRA_IDS)
            shell:
                "mkdir -p input_reads && "
                "ln -sf $(readlink -f {input.r1}) {output.r1} && "
                "ln -sf $(readlink -f {input.r2}) {output.r2}"
    
    rule trimming_reads:
        input:
            r1 = "input_reads/{sample}_1.fastq",
            r2 = "input_reads/{sample}_2.fastq"
        output:
            p1 = "clean_data/{sample}_1P.fastq",
            p2 = "clean_data/{sample}_2P.fastq"
        shell:
            """
            mkdir -p clean_data && 
            trimmomatic PE -threads {threads} {input.r1} {input.r2} \
            {output.p1} clean_data/{wildcards.sample}_1U.fastq \
            {output.p2} clean_data/{wildcards.sample}_2U.fastq \
            LEADING:20 TRAILING:20 MINLEN:20 \
            """

    rule genome_assembly_sra:
        input:
            r1 = "clean_data/{sample}_1P.fastq",
            r2 = "clean_data/{sample}_2P.fastq"
        output:
            "assembly/{sample}/contigs.fasta"
        wildcard_constraints:
            sample="|".join(SRA_IDS)
        threads: 8
        shell:
            "mkdir -p $(dirname {output}) && "
            "spades.py -t {threads} -1 {input.r1} -2 {input.r2} -o $(dirname {output})"

if ASSEMBLY_IDS:
    if ASSEMBLY_SOURCE == "online":
        rule download_assembly:
            output:
                "assembly/{sample}/contigs.fasta"
            wildcard_constraints:
                sample="|".join(ASSEMBLY_IDS)
            shell:
                "mkdir -p $(dirname {output}) && "
                "datasets download genome accession {wildcards.sample} --filename {wildcards.sample}.zip && "
                "unzip -o {wildcards.sample}.zip -d $(dirname {output}) && "
                "find $(dirname {output}) -name '*.fna' -exec mv {{}} {output} \; && "
                "rm {wildcards.sample}.zip"
        
    elif ASSEMBLY_SOURCE == "offline":
        rule local_assembly_input:
            input:
                lambda wildcards: os.path.join(ASSEMBLY_LOCAL_DIR, f"{wildcards.sample}.fasta")
            output:
                "assembly/{sample}/contigs.fasta"
            wildcard_constraints:
                sample="|".join(ASSEMBLY_IDS)
            shell:
                "mkdir -p $(dirname {output}) && "
                "cp {input} {output}"
rule run_quast:
    input:
        "assembly/{sample}/contigs.fasta"
    output:
        "quast_output/{sample}/report.tsv"
    shell:
        "mkdir -p $(dirname {output}) && "
        "quast {input} -o $(dirname {output})"

rule run_prokka:
    input:
        "assembly/{sample}/contigs.fasta"
    output:
        gff = "prokka_output/{sample}/{sample}.gff",
        tsv = "prokka_output/{sample}/{sample}.tsv",
    shell:
        "mkdir -p $(dirname {output.gff}) && "
        "prokka --force --outdir $(dirname {output.gff}) --prefix {wildcards.sample} {input}"

rule amine_search:
    input:
        "prokka_output/{sample}/{sample}.tsv"
    output:
        "decarboxylases/amines_{sample}/amines_{sample}.tsv"
    shell:
        """
        mkdir -p $(dirname {output}) && 
        awk -F '\t' 'tolower($4) ~ /(gcv|gdc|aldc|tdca?b?|tdc[-]?1|ddc|aadc|tydc|tyrdc|mfn[ap]|adc|adia?|hdc|trp[ab]|cyp71p1|odc|spe[acdefb]|ldcc?|cad[ab]|decarboxylase)/ {{
            print $1 "\t" $4 "\t" $7
        }}' {input} > {output}
        """

rule presence_absence_matrix:
    input:
        sample_files=expand("decarboxylases/amines_{sample}/amines_{sample}.tsv", sample=ALL_IDS)
    output:
        "summary/presence_absence_matrix.tsv"
    run:
        import os
        import pandas as pd
        import re
        from collections import defaultdict
        
        os.makedirs("summary", exist_ok=True)
        gene_pattern = re.compile(
            r'(gcv|gdc|aldc|tdca?b?|tdc[-]?1|ddc|aadc|tydc|tyrdc|mfn[ap]|adc|adia?|hdc|trp[ab]|cyp71p1|odc|spe[acdefb]|ldcc?|cad[ab]|decarboxylase)',
            re.IGNORECASE
        )
        
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

def get_input_type(wildcards):
    if wildcards.sample in ASSEMBLY_IDS:
        return "assembly"
    elif wildcards.sample in SRA_IDS:
        return "reads"
    else:
        raise ValueError(f"Unknown sample type: {wildcards.sample}")

rule generate_sample_reports:
    input:
        decarboxylases="decarboxylases/amines_{sample}/amines_{sample}.tsv",
        quast_report="quast_output/{sample}/report.tsv",
        assembly="assembly/{sample}/contigs.fasta"
    output:
        "reports/sample_{sample}.html"
    params:
        input_type=lambda wildcards: get_input_type(wildcards)
    run:
        import pandas as pd
        import os
        
        sample = wildcards.sample
        os.makedirs("reports", exist_ok=True)
        
        with open(output[0], 'w') as f:
            # HTML-head
            f.write(f"""
            <!DOCTYPE html>
            <html>
            <head>
                <meta charset="UTF-8">
                <title>Report for {sample}</title>
                <style>
                    body {{ font-family: Arial, sans-serif; margin: 20px; }}
                    h1, h2 {{ color: #2c3e50; }}
                    table {{ border-collapse: collapse; margin: 15px 0; }}
                    th, td {{ border: 1px solid #ddd; padding: 6px; text-align: left; }}
                    th {{ background-color: #f2f2f2; }}
                    .no-genes {{ color: #4CAF50; font-weight: bold; }}
                    .warning {{ color: #FF9800; }}
                </style>
            </head>
            <body>
            <h1>Report for Sample: {sample}</h1>
            """)
            
            # genes of interest
            f.write(f"<h2>Genes of Interest</h2><p>{', '.join(GENES)}</p>")
            
            # detected genes
            f.write("<h2>Detected Decarboxylase Genes</h2>")
            if os.path.exists(input.decarboxylases) and os.path.getsize(input.decarboxylases) > 0:
                try:
                    amines_df = pd.read_csv(input.decarboxylases, sep='\t', header=None, names=["locus_tag", "gene", "product"])
                    f.write(amines_df.to_html(index=False, header=False))
                except Exception as e:
                    f.write(f'<p class="warning">Error reading amines file: {str(e)}</p>')
            else:
                f.write('<p class="no-genes">🟢 No biogenic amine-related genes detected</p>')
            
            # QUAST report - only show for SRA samples (reads)
            if params.input_type == "reads":
                f.write("<h2>Assembly Quality Metrics (QUAST)</h2>")
                if os.path.exists(input.quast_report) and os.path.getsize(input.quast_report) > 0:
                    try:
                        quast_df = pd.read_csv(input.quast_report, sep='\t')
                        f.write(quast_df.to_html(index=False))
                    except Exception as e:
                        f.write(f'<p class="warning">Error reading QUAST report: {str(e)}</p>')
                else:
                    f.write('<p class="warning">QUAST report not available</p>')
            
            # assembly info
            f.write("<h2>Assembly Information</h2>")
            if os.path.exists(input.assembly):
                try:
                    with open(input.assembly, 'r') as fasta:
                        contig_count = sum(1 for line in fasta if line.startswith('>'))
                    f.write(f"<p>Contig count: {contig_count}</p>")
                except Exception as e:
                    f.write(f'<p class="warning">Error reading assembly: {str(e)}</p>')
            else:
                f.write('<p class="warning">Assembly file not available</p>')
            
            f.write("</body></html>")

rule aggregate_reports:
    input:
        sample_reports = expand("reports/sample_{sample}.html", sample=ALL_IDS),
        heatmap_png = "summary/presence_absence_heatmap.png",
        matrix_tsv = "summary/presence_absence_matrix.tsv"
    output:
        "reports/full_report.html"
    run:
        import os
        import shutil
        from datetime import datetime
        
        os.makedirs("reports", exist_ok=True)
        
        with open(output[0], 'w') as out_f:
            out_f.write("""
            <!DOCTYPE html>
            <html>
            <head>
                <meta charset="UTF-8">
                <title>Biogenic Amine Analysis Report</title>
                <style>
                    body { font-family: Arial, sans-serif; margin: 20px; }
                    h1, h2 { color: #2c3e50; }
                    .sample-report { margin: 25px 0; padding: 15px; border: 1px solid #eee; }
                    .heatmap-container { text-align: center; margin: 30px 0; }
                    .matrix-link { margin: 20px 0; }
                    .tools-table { margin: 30px 0; width: 100%; }
                    img { max-width: 80%; }
                </style>
            </head>
            <body>
            <h1>Biogenic Amine Analysis Report</h1>
            """)
            
            out_f.write(f"<p><strong>Date:</strong> {datetime.now().strftime('%Y-%m-%d')}</p>")
            out_f.write(f"<p><strong>Samples analyzed:</strong> {len(ALL_IDS)}</p>")
            
            # heatmap
            out_f.write("""
            <div class='heatmap-container'>
                <h2>Summary: Presence/Absence Heatmap</h2>
                <img src='heatmap.png' alt='Gene Presence Heatmap'>
            </div>
            """)
            
            # matrix
            out_f.write(f"""
            <div class='matrix-link'>
                <h2>Gene Presence/Absence Matrix</h2>
                <p>Download full table: <a href='summary/presence_absence_matrix.tsv'>matrix.tsv</a></p>
            </div>
            """)
            
            # sample reports
            out_f.write("<h2>Sample Reports</h2>")
            for report_path in input.sample_reports:
                sample = os.path.basename(report_path).replace('sample_', '').replace('.html', '')
                with open(report_path, 'r') as in_f:
                    content = in_f.read()
                    start = content.find('<body>') + 6
                    end = content.find('</body>')
                    if start > 5 and end > start:
                        sample_content = content[start:end].strip()
                        out_f.write(f"<div class='sample-report'><h3>Sample: {sample}</h3>{sample_content}</div>")
            
            # tools
            out_f.write("""
            <h2>Pipeline Tools</h2>
            <table class='tools-table'>
                <tr><th>Tool</th><th>Version</th><th>Purpose</th></tr>
                <tr><td>SPAdes</td><td>3.15.5</td><td>Genome assembly</td></tr>
                <tr><td>QUAST</td><td>5.2.0</td><td>Quality assessment</td></tr>
                <tr><td>Prokka</td><td>1.14.6</td><td>Annotation</td></tr>
                <tr><td>Trimmomatic</td><td>0.39</td><td>Read trimming</td></tr>
            </table>
            </body></html>
            """)
        
        shutil.copy(input.heatmap_png, "reports/heatmap.png")
        shutil.copy(input.matrix_tsv, "reports/matrix.tsv")