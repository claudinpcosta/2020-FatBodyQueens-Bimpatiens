##################################################################### 
#Read mapping
# Analysis by Joshua P. Der
# <jder@fullerton.edu>
# Last updated: 18 July 2016
# Data correspond to Costa et al. 2020 
#####################################################################


#==================
## Raw Illumina data is available from the NCBI Sequence Read Archive under the accession number (...). 


#==================
## FASTQC was run on the forward and reverse reads in two runs to 

module load genomics/all

cd ./160614_K00180_0207_AHF277BBXX/Reports/FASTQC

fastqc -o ./160614_K00180_0207_AHF277BBXX/Reports/FASTQC/ --noextract -t 8 ./160614_K00180_0207_AHF277BBXX/raw_data/Q*R1_001.fastq.gz

fastqc -o ./160614_K00180_0207_AHF277BBXX/Reports/FASTQC/ --noextract -t 8 ./160614_K00180_0207_AHF277BBXX/raw_data/Q*R2_001.fastq.gz


#===================
## Run Trimmomatic:

module load genomics/trimmomatic
# adapters:
# https://github.com/timflutre/trimmomatic/blob/master/adapters/TruSeq3-PE-2.fa

cd ./160614_K00180_0207_AHF277BBXX/raw_data
for i in Q*_R1_001.fastq.gz
do
    sample="$(basename $i _R1_001.fastq.gz)"
    echo working with $sample
    java -jar ./Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 -threads 8 ${sample}_R1_001.fastq.gz ${sample}_R2_001.fastq.gz ../trim+filter/${sample}.trim.R1.fq.gz ../trim+filter/${sample}.trim.U1.fq.gz ../trim+filter/${sample}.trim.R2.fq.gz ../trim+filter/${sample}.trim.U2.fq.gz ILLUMINACLIP:./Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 
done

cd ../trim+filter
for i in *.U1.fq.gz
do
	sample="$(basename $i .U1.fq.gz)"
	cat ${sample}.U1.fq.gz ${sample}.U2.fq.gz > ${sample}.SE.fq.gz
	rm ${sample}.U1.fq.gz ${sample}.U2.fq.gz
done


#=========================
## Combine files by sample ID (i.e. merge the three lanes):


cd ./160614_K00180_0207_AHF277BBXX/trim+filter

for i in Q*_L003.trim.R1.fq.gz
do
    sample="$(basename $i _L003.trim.R1.fq.gz)"
    cat ${sample}_L003.trim.R1.fq.gz ${sample}_L004.trim.R1.fq.gz ${sample}_L005.trim.R1.fq.gz > ${sample}.trim.R1.fq.gz &
    cat ${sample}_L003.trim.R2.fq.gz ${sample}_L004.trim.R2.fq.gz ${sample}_L005.trim.R2.fq.gz > ${sample}.trim.R2.fq.gz &
    cat ${sample}_L003.trim.SE.fq.gz ${sample}_L004.trim.SE.fq.gz ${sample}_L005.trim.SE.fq.gz > ${sample}.trim.SE.fq.gz &
done

###ARGH files got corrupted during cat
#checking file lengths:

for i in *R1.fq.gz
do
x="$(basename $i R1.fq.gz)"
echo $x
zcat ${x}R1.fq.gz |wc -l
zcat ${x}R2.fq.gz |wc -l
echo
done

# counting reads: 
for i in *fq
do
lines="$(cat $i | wc -l)"
reads=$((lines/4))
printf "$i\t$reads\n"
done



#### Recombining files
for i in *_L003.trim.R1.fq.gz
do
    sample="$(basename $i _L003.trim.R1.fq.gz)"
    zcat ${sample}_L003.trim.R1.fq.gz ${sample}_L004.trim.R1.fq.gz ${sample}_L005.trim.R1.fq.gz > ${sample}.trim.R1.fq
    zcat ${sample}_L003.trim.R2.fq.gz ${sample}_L004.trim.R2.fq.gz ${sample}_L005.trim.R2.fq.gz > ${sample}.trim.R2.fq
    zcat ${sample}_L003.trim.SE.fq.gz ${sample}_L004.trim.SE.fq.gz ${sample}_L005.trim.SE.fq.gz > ${sample}.trim.SE.fq
done





#==========================================
## Build the HiSat2 index

module load genomics/hisat2
module load python/2.7.8

cd ./160614_K00180_0207_AHF277BBXX/Bimp_2.0_Sadd/Bimp_2.0.hisat2/

hisat2-build -p 32 --ss ../bimp_OGSv1.0.ss --exon ../bimp_OGSv1.0.exons -f ../Bimp_2.0_scaffolds.fa Bimp_2.0.h2i

echo "DONE"
date


#==========================================
### Map reads using hisat

cd ./160614_K00180_0207_AHF277BBXX/mapping
module load genomics/hisat2

for i in *_L003.trim.R1.fq.gz
do
    sample="$(basename $i _L003.trim.R1.fq.gz)"
    echo $sample
    hisat2 -p 40 -q -x ../Bimp_2.0_Sadd/Bimp_2.0.ht2/Bimp_2.0.hti -1 ${sample}_L003.trim.R1.fq.gz,${sample}_L004.trim.R1.fq.gz,${sample}_L005.trim.R1.fq.gz -2 ${sample}_L003.trim.R2.fq.gz,${sample}_L004.trim.R2.fq.gz,${sample}_L005.trim.R2.fq.gz -U ${sample}_L003.trim.SE.fq.gz,${sample}_L004.trim.SE.fq.gz,${sample}_L005.trim.SE.fq.gz -S ${sample}.sam 2>${sample}.hisat2.err
done

# Extract alignment rates:
perl -e 'foreach $file (@ARGV){ open (FILE, $file); while(<FILE>){if(/(\d+\.\d+)% overall alignment rate/){print "$file\t$1\n"} } } ' *.err > overall_alignment_rate.txt


### sam to bam conversion
# I initially ran this to run stringtie for abundances, but stringtie only outputs normalized 
# expression levels (fpkm or transcripts per million) and I need counts to use DESeq and/or EdgeR. 
cd ./160614_K00180_0207_AHF277BBXX/mapping
module load genomics/samtools
for i in *.sam
do
    sample="$(basename $i .sam)"
    echo $sample
    samtools view -@ 40 -u ${sample}.sam | samtools sort -@ 40 -o ${sample}.bam -
done

###==============================================
# Run featureCounts to identify extract and compile read counts found in gene features that map to exons.

cd ./160614_K00180_0207_AHF277BBXX/quantitation
# run featureCounts
./subread-1.5.0-p3-Linux-x86_64/bin/featureCounts -g gene_id -t exon -T 32 -p -a ../Bimp_2.0_Sadd/bimp_OGSv1.0.gtf -o gene_counts.txt ./*.bam



####=======================================
# Preparing files for DE analysis
# I modified the header line in gene_counts.txt file to clean up the sample names (I dropped things like _S##.bam).
# Then I extracted experimental design from the sample names using the following perl onliner:

perl -e 'print "\tqueen\tcolony\tage\tdiet\n"; while(<>){chomp; $id = $_; $id =~ /(Q\d+)_(R\d+)_(\d+)d(\d*)_S/; if($4 eq ""){print "$id\t$1\t$2\t$3\t0\n";} else{print "$id\t$1\t$2\t$3\t$4\n";}}' samples.txt > design.txt 

# The above model works, but it actually groups the 0 day queens with the 0% diet treatment, which is not correct. When I code the 0 day bees as having a different diet treatment (NA = missing data, or "none" being distinct from the other treatments, the model is not fully specified (age fully predicts the diet condition for 0 day bees). This results in model specification failing and I need to consult a statistician about how to specify the experimental design.

# The experimental design needs to account for no diet treatment in the 0 day queens. 
# The following oneliner combines the age+diet combination into a single factor. I'm not sure if this is best stats design, but it produced a full rank model.
perl -e 'while(<>){chomp; @a=split /\t/; print "$a[0]\t$a[2]\t$a[3]". "_". "$a[4]\n"}' design.v3.txt > design.v4.txt 


###======================
# DEG analysis, see R coe (...). 



