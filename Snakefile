# used fasterq-dump./ # with # being each SRA number 
#datasets download genome accession GCF_000845245.1 --include gff3,rna,cds,protein,genome,seq-report

#used this command line to download the cds from NCBI
#datasets download genome accession GCF_000845245.1 --include gff3,rna,cds,protein,genome,seq-report

#### Step 2 ####

rule downloading_genome: 
    output:
        directory("dataset/GCF_000845245.1"),

    shell:
        """
        mkdir -p dataset/GCF_000845245.1
        datasets download genome accession GCF_000845245.1 --include cds,protein,genome --filename "dataset/ncbi_dataset.zip"
        unzip -o "dataset/ncbi_dataset.zip" -d dataset/
        cp dataset/ncbi_dataset/data/GCF_000845245.1/cds_from_genomic.fna {output}
        cp dataset/ncbi_dataset/data/GCF_000845245.1/GCF_000845245.1_ViralProj14559_genomic.fna {output}
        cp dataset/ncbi_dataset/data/GCF_000845245.1/protein.faa {output}
    
        """
rule make_list: 
    input:
        "dataset/GCF_000845245.1/cds_from_genomic.fna"
    output: 
        "dataset/output.fasta",
        "dataset/PipelineReport.txt"
    run:
        from Bio import SeqIO
        seq = []
        id_list = []
        x = []
        seq_pair= []
        protein_list = []

        for record in SeqIO.parse(input[0], "fasta"):
            id_list.append(record.id)
            seq.append(str(record.seq))

        #print(seq)
        for i in id_list: 
            x = list(i)
            #print(x)
            for y in x:
                if y == "Y":
                    pos = x.index(y)
                    protein_name = (x[pos:])
                    protein_name = ''.join(protein_name)
                    protein_name = "".join([">", protein_name])
                    protein_list.append(protein_name)
        count = 0
        for w in range(len(seq)):
            seq_pairs = "".join([protein_list[w] +"\n", seq[w]])
            seq_pair.append(seq_pairs)
            count += 1 

        #print(count)

        with open (output[0], "w") as outfile:
            outfile.write("\n".join(seq_pair))

        with open (output[1], "w") as outfile:
            outfile.write(str(f"The HCMV genome (GCF_000845245.1) has {count} CDS."))

rule kallisto:
    input: 
        "dataset/output.fasta"
    output:
        "dataset/kallisto_output.idx",
    shell:
       "kallisto index -i {output} {input}"

#### STEP 3 ####
rule kallisto_quantify:
    input: 
        idx = "dataset/kallisto_output.idx",
        d1_2in_f = "SRR5660030_1.fastq",
        d1_2in_r = "SRR5660030_2.fastq",
        d1_6in_f = "SRR5660033_1.fastq",
        d1_6in_r = "SRR5660033_2.fastq",
        d2_2in_f = "SRR5660044_1.fastq",
        d2_2in_r = "SRR5660044_2.fastq",
        d2_6in_f = "SRR5660045_1.fastq",
        d2_6in_r = "SRR5660045_2.fastq"

    output:
        d1_2 = directory("results/SRR5660030"),
        d1_6 = directory("results/SRR5660033"),
        d2_2 = directory("results/SRR5660044"),
        d2_6 = directory("results/SRR5660045")

    shell:
        """
        kallisto quant -i {input.idx} -o {output.d1_2} -b 30 {input.d1_2in_f} {input.d1_2in_r}
        kallisto quant -i {input.idx} -o {output.d1_6} -b 30 {input.d1_6in_f} {input.d1_6in_r}
        kallisto quant -i {input.idx} -o {output.d2_2} -b 30 {input.d2_2in_f} {input.d2_2in_r}
        kallisto quant -i {input.idx} -o {output.d2_6} -b 30 {input.d2_6in_f} {input.d2_6in_r}
        """
rule make_sleuth_table: 
    input: 
        d1_2 = "results/SRR5660030/abundance.h5",
        d1_6 = "results/SRR5660033/abundance.h5",
        d2_2 = "results/SRR5660044/abundance.h5",
        d2_6 = "results/SRR5660045/abundance.h5"
    output: 
       table = "results/sleuth.tsv"
    run: 
        import pandas as pd
        #time names (horizontal)
        matrix_frame = {
        "sample":["SRR5660030","SRR5660033","SRR5660044", "SRR5660045"],
        "condition":["2dpi", "6dpi", "2dpi", "6dpi"],
        "path":[input.d1_2, input.d1_6, input.d2_2, input.d2_6]
        }

        sleuth_matrix = pd.DataFrame(matrix_frame)

        with open (output.table, "w") as outfile:
            outfile.write(str(sleuth_matrix))

rule sleuth_run: 
    input: 
        table = "results/sleuth.tsv"
    output: 
        sleuth_out = "results/sleuth_output.txt"
    script:
        "sleuth.R"

rule add_to_pipeline_txt:
    input: 
        "results/sleuth_output.txt"
    output:
        "dataset/PipelineReport.txt"
    shell:
        "cat {input} >> {output}"


#### step 4 ####
sample = ["SRR5660030","SRR5660033","SRR5660044","SRR5660045"]

rule all:
    input: 
        expand("results/bowtie_data_out/{sample}_map.sam", sample = sample)
    
rule bowtie2_genome_index: 
    input: 
        genome = "dataset/ncbi_dataset/data/GCF_000845245.1/GCF_000845245.1_ViralProj14559_genomic.fna",
    
    output: 
        directory("results/bowtie_data/GCF_index")
    shell: 
       """
       mkdir -p {output} &&
       bowtie2-build {input.genome} {output}/GCF_index
       """

rule running_bowtie_2:
    input:
        genome_results = "results/bowtie_data/GCF_index/GCF_index",
        forward_reads = "{sample}_1.fastq",
        reverse_reads = "{sample}_2.fastq"
    output: 

        "results/bowtie_data_out/{sample}.sam"
    shell: 
        """
        mkdir -p results/bowtie_data_out &&
        bowtie2 --quiet -x {input.genome_results} -1 {input.forward_reads} -2 {input.reverse_reads} -S {output}
        """