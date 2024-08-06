#!/bin/bash -l

############# SLURM SETTINGS #############
#SBATCH --account=project0026 
#SBATCH --job-name=P1 
#SBATCH --partition=nodes 
#SBATCH --time=0-04:00:00 
#SBATCH --mem=4G 
#SBATCH --nodes=1 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=2 
#SBATCH --ntasks-per-node=1 
#SBATCH --mail-user=2938235s@student.gla.ac.uk 
#SBATCH --mail-type=END  
#SBATCH --mail-type=FAIL
############# CODE #############

data=/mnt/autofs/data/userdata/project0026/pathogenPolyomicsData/DNAseq/
cp -r /mnt/autofs/data/userdata/project0026/pathogenPolyomicsData/Reference .
referencedir=Reference/

outputdir=P1output/
mkdir -p ${outputdir}

fastqcdir=${outputdir}fastqc/
mkdir -p ${fastqcdir}

trimfastqcdir=${outputdir}fastqc/trimmed
mkdir -p ${trimfastqcdir}

trimdir=${outputdir}trimgalore/
mkdir -p ${trimdir}

indexdir=${outputdir}index/
mkdir -p ${indexdir}

bowtiedir=${outputdir}bowtie/
mkdir -p ${bowtiedir}

bamdir=${outputdir}bam/
mkdir -p ${bamdir}

vcfdir=${outputdir}vcf/
mkdir -p ${vcfdir}

for sample in AmpB_1 AmpB_2 WT_1 WT_2 
do
	sample=Lmex${sample}
	fastq=${data}${sample}.fastq.gz
	fastqc -o ${fastqcdir} ${fastq}
done

bowtie2-build ${referencedir}TriTrypDB-25_LmexicanaMHOMGT2001U1103.fa ${indexdir}Lmex

for sample in AmpB WT
do
	sample=Lmex${sample}
	trim_galore --paired --phred64 -o ${trimdir} ${data}${sample}_1.fastq.gz ${data}${sample}_2.fastq.gz
	bowtie2 --phred64 -q -S ${bowtiedir}${sample}.sam -x ${indexdir}Lmex -1 ${trimdir}${sample}_1_val_1.fq.gz -2 ${trimdir}${sample}_2_val_2.fq.gz 2> ${outputdir}${sample}_bowtie.log
	samtools view -b -u ${bowtiedir}${sample}.sam |	samtools sort -o ${bamdir}${sample}_sorted.bam
	samtools index ${bamdir}${sample}_sorted.bam
	bamCoverage -b ${bamdir}${sample}_sorted.bam -o ${bamdir}${sample}_sorted.bam.bw
	bamCoverage -b ${bamdir}${sample}_sorted.bam -of bedgraph -o ${bamdir}${sample}_sorted.bam.bedgraph
	samtools flagstat ${bamdir}${sample}_sorted.bam
done

for sample in AmpB_1_val_1 AmpB_2_val_2 WT_1_val_1 WT_2_val_2 
do
	sample=Lmex${sample}
	fastq=${trimdir}${sample}.fq.gz
	fastqc -o ${trimfastqcdir} ${fastq}
done

bamaddrg -b ${bamdir}LmexWT_sorted.bam -s WT -b ${bamdir}LmexAmpB_sorted.bam -s AmpB | freebayes --stdin -p 2 -f ${referencedir}TriTrypDB-25_LmexicanaMHOMGT2001U1103.fa >${vcfdir}Lmex.vcf

vcffilter -f "QUAL > 20" -f "TYPE = snp" ${vcfdir}Lmex.vcf >${vcfdir}Lmex_qualfilt.vcf

for sample in WT AmpB
do
	bcftools view -c 1 -Oz -s ${sample} -o ${vcfdir}Lmex_${sample}_filt.vcf.gz ${vcfdir}Lmex_qualfilt.vcf
	bcftools index ${vcfdir}Lmex_${sample}_filt.vcf.gz
done

isecdir=${vcfdir}bcfisec/
mkdir -p ${isecdir}

bcftools isec -p ${isecdir} -Ov ${vcfdir}Lmex_WT_filt.vcf.gz ${vcfdir}Lmex_AmpB_filt.vcf.gz
mv ${isecdir}0000.vcf ${isecdir}WT_only.vcf
mv ${isecdir}0001.vcf ${isecdir}AmpB_only.vcf
mv ${isecdir}0002.vcf ${isecdir}WT_int_AmpB.vcf
mv ${isecdir}0003.vcf ${isecdir}AmpB_int_WT.vcf

snpdir=${vcfdir}snp/
mkdir -p ${snpdir}

#snpEff build -c SnpEff.config -gff3 -noCheckCds -noCheckProtein -v Lmex
snpEff -Xmx4g -c SnpEff.config -no-intron -no SPLICE_SITE_REGION Lmex ${isecdir}AmpB_only.vcf > ${snpdir}AmpB_only_anno.vcf
cat ${snpdir}AmpB_only_anno.vcf | vcfEffOnePerLine.pl | SnpSift extractFields - CHROM POS REF ALT "ANN[*].IMPACT" "ANN[*].EFFECT" "ANN[*].GENE" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "GEN[*].GT" | grep 'HIGH\|MODERATE' | grep '1/1' > ${snpdir}AmpB_sifted.tsv

