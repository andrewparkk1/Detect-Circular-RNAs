#!/bin/bash
############################################################
# Help                                                     #
############################################################
Help() {
   # Display Help
   echo "Syntax: scriptTemplate [-g|h|v]"
   echo "options:"
   echo "g     Print the GPL license notification."
   echo "h     Print this Help."
   echo "v     Print software version and exit."
   echo
}

############################################################
############################################################
# Main program                                             #
############################################################
############################################################

############################################################
# Process the input options. Add options as needed.        #
############################################################
# Get the options
while :; do
    case $1 in
        -h|-\?|--help)
            Help    # Display a usage synopsis.
            exit;;
        -rnas1|-1) 
            shift
            i1=$1;;
        -rnas2|-2) 
            shift
            i2=$1;;
        -ribos1|-3)
            shift
            i3=$1;;
        -ribos2|-4)
            shift 
            i4=$1;;
        -trna)
            shift
            trna=$1;;
        -rrna)
            shift
            rrna=$1;;
        -fa)
            shift
            fa=$1;;
        -gtf)
            shift
            gtf=$1;;
        -?*)
            printf 'WARN: Unknown option (ignored): %s\n' "$1" >&2
            ;;
        *)               # Default case: No more options, so break out of the loop.
            break
    esac
    shift
done


############################################################
#                               EXTRA                      #
############################################################
# fa="Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
# gtf="Homo_sapiens.GRCh38.105.chr.gtf"
# i1="SRR12354633 SRR12354634"``

subs=$(ls -d */)
mkdir data
mkdir data/folders
mkdir data/sra
mkdir data/rna-seq
mkdir data/rna-seq/sample1
mkdir data/rna-seq/sample2
mkdir data/ribo-seq
mkdir genome-index
mkdir fastqc
mkdir fastqc/rna
mkdir fastqc/ribo


# i1="SRR12354633 SRR12354634 SRR12354635"
# i2="SRR12354636 SRR12354637 SRR12354638"
# i3="SRR12354645 SRR12354646 SRR12354647"
# i4="SRR12354649 SRR12354650 SRR12354651"

r1=($i1)
r2=($i2)
r3=($i3)
r4=($i4)

# #move SRR folders to separate folder
# #move SRR files to separate folder
for i in $subs; do
  mv $PWD/$i $PWD/data/folders
  mv $PWD/data/folders/$i/* $PWD/data/sra
done
rm -r data/folders


# #loop through each set 
# #in each set, fastq-dump, rename, fastqc


j=1
for i in "${r1[@]}"
do
    # fastq-dump $PWD/data/sra/$i.sra --outdir $PWD/data/rna-seq
    # mv $PWD/data/rna-seq/$i.fastq $PWD/data/rna-seq/sample1/rna_sample1_rep$j.fastq
    fastqc $PWD/data/rna-seq/sample1/rna_sample1_rep$j.fastq --outdir=$PWD/fastqc/rna
    ((j=j+1))
done

j=1
for i in "${r2[@]}"
do
    # fastq-dump $PWD/data/sra/$i.sra --outdir $PWD/data/rna-seq
    # mv $PWD/data/rna-seq/$i.fastq $PWD/data/rna-seq/sample2/rna_sample2_rep$j.fastq
    fastqc $PWD/data/rna-seq/sample2/rna_sample2_rep$j.fastq --outdir=$PWD/fastqc/rna
    ((j=j+1))
done

j=1
for i in "${r3[@]}"
do
    # fastq-dump $PWD/data/sra/$i.sra --outdir $PWD/data/ribo-seq
    # mv $PWD/data/ribo-seq/$i.fastq $PWD/data/ribo-seq/ribo_sample1_rep$j.fastq
    fastqc $PWD/data/ribo-seq/ribo_sample1_rep$j.fastq --outdir=$PWD/fastqc/ribo
    ((j=j+1))
done

j=1
for i in "${r4[@]}"
do
    # fastq-dump $PWD/data/sra/$i.sra --outdir $PWD/data/ribo-seq
    # mv $PWD/data/ribo-seq/$i.fastq $PWD/data/ribo-seq/ribo_sample2_rep$j.fastq
    fastqc $PWD/data/ribo-seq/ribo_sample2_rep$j.fastq --outdir=$PWD/fastqc/ribo
    ((j=j+1))
done

multiqc $PWD/fastqc/rna -o $PWD/fastqc/rna
multiqc $PWD/fastqc/ribo -o $PWD/fastqc/ribo

# Download Ensemble Mouse/Human FASTA (.fa) and GTF (.gtf)
# wget http://ftp.ensembl.org/pub/release-104/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
# wget http://ftp.ensembl.org/pub/release-104/gtf/mus_musculus/Mus_musculus.GRCm39.104.chr.gtf.gz

# Create genome index
# http://useast.ensembl.org/Mus_musculus/Info/Index
STAR --runThreadN 32 --runMode genomeGenerate --genomeDir genome-index --genomeFastaFiles $fa --sjdbGTFfile $gtf --sjdbOverhang 100 --genomeChrBinNbits 15

# Detect circRNAs
rnaS1=$(ls -d $PWD/data/rna-seq/sample1/* | xargs echo | sed 's/ /,/g')
rnaS2=$(ls -d $PWD/data/rna-seq/sample2/* | xargs echo | sed 's/ /,/g')

# $cd seekCRIT
# $(pip3 install -r Prerequisites.txt)
# $(python setup.py install)
mkdir seekCRIT-output
python3 $PWD/seekCRIT/seekCRIT/seekCRIT.py -o seekCRIT-output -t SE -fa $fa -gtf $gtf --threadNumber 32 --genomeIndex genome-index -s1 $rnaS1 -s2 $rnaS2 -ref refseq.txt

# python3 seekCRIT.py -o seekCRIT-output -t PE -fa $fa -gtf $gtf --threadNumber 32 --genomeIndex genome-index 
# -s1 wt1_R1.fastq:wt1_R2.fastq,wt2_R1.fastq:wt2_R2.fastq,wt3_R1.fastq:wt3_R2.fastq 
# -s2 mettl1_R1.fastq:mettl1_R2.fastq,mettl2_R1.fastq:mettl2_R2.fastq,mettl3_R1.fastq:mettl3_R2.fastq 

# Create bed file to get sequence 
$(awk -F "\t" 'OFS="\t" {print $1,$2,$3,$10,"0",$4}' circRNAs.pVal.FDR.txt > circRNA.bed)
$(sed '1d' circRNA.bed > c.bed)

# # $(brew tap homebrew/science)
# # $(brew install bedtools)

# # # # Use bedtools to retrieve exon sequence of circRNA
$(bedtools getfasta -s -fi $fa -bed c.bed > sequence.txt)

# # # # Create 100bp-sequence per circRNA using sequence.txt.
python bsj-seq.py

# # Create index for circRNAs using the seekCRIT output
$(mkdir circ-genome-index)
STAR --runThreadN 32 --runMode genomeGenerate --genomeDir circ-genome-index --genomeFastaFiles 100bp-circRNA.fa

# Obtain tRNA and rRNA data
mkdir tRNA-index
mv $trna tRNA-index
mkdir rRNA-index
mv $rrna rRNA-index
STAR --runThreadN 32 --runMode genomeGenerate --genomeDir trna-index --genomeFastaFiles tRNA-index/$trna --genomeSAindexNbases 6
STAR --runThreadN 32 --runMode genomeGenerate --genomeDir rrna-index --genomeFastaFiles rRNA-index/$rrna --genomeSAindexNbases 6

# Trim Ribo-seq of its adaptor sequences using Trim Galore
# May need to install cut_adaptor separately
# leave everything less than 40bp
mkdir $dir/data/trimmed-ribo-seq
trim_galore -o $PWD/$dir/data/trimmed-ribo-seq --max_length 40 # $PWD/data/ribo-seq/sample1/ribo1_sample1.fastq 

# Run FastQC on TRIMMED Ribo-seq
fastqc $PWD/data/trimmed-ribo-seq/rna1_sample2_trimmed.fq --outdir=$PWD/fastqc/rna

# # Map riboseq to tRNA, rRNA, chr, circular
# # Make sure the unmapped reads have separate output
# # Map ribo to trna genome index
STAR --runThreadN 32 --genomeDir $PWD/tRNA-index --outSAMunmapped Within --outReadsUnmapped Fastx --outFileNamePrefix tRNA-wt1- --readFilesIn $PWD/data/trimmed-ribo-seq/*

# # Map unmapped-trna-genome to rrna genome index
STAR --runThreadN 32 --genomeDir rRNA-genome-index --alignEndsType EndToEnd --outSAMunmapped Within --outReadsUnmapped Fastx --outFileNamePrefix tRNA- --readFilesIn tRNA-wt1-unmapped.fa

# # Map unmapped-rrna-trna-genome to mouse genome index
STAR --runThreadN 32 --genomeDir mouse-genome-index --alignEndsType EndToEnd --outSAMunmapped Within --outReadsUnmapped Fastx --outFileNamePrefix tRNA- --readFilesIn <unmapped rRNA fastq>

# # Map unmapped-chr-rrna-trna-genome to circularRNA index
STAR --runThreadN 32 --genomeDir circular-genome-index --alignEndsType EndToEnd --outFileNamePrefix circular --readFilesIn <unmapped chr fastq>
