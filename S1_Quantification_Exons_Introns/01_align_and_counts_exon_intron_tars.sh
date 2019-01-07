#!/bin/bash
#sbatch --array 1-970  --mem=20000 -J ALN_CNT -o "/pasteur/projets/policy01/evo_immuno_pop/Maxime/Splicing/HISAT2/JobOutput/JobOutput_%a.log" /pasteur/homes/mrotival/03_Analysis/Splicing/Papier_Splicing/ReviewNcomm/all_scripts_to_share/align_and_counts_exon_intron_tars.sh



USER='mrotival'
SCRIPTS="/pasteur/projets/policy01/evo_immuno_pop/Maxime/Splicing/HISAT2/hisat2-master/"
HISATDIR="/pasteur/projets/policy01/evo_immuno_pop/Maxime/Splicing/HISAT2/"
TMPDIR="/pasteur/scratch/${USER}"

SAMPLE_NUM=${SLURM_ARRAY_TASK_ID}
SAMPLEDIR=${HISATDIR}

source /local/gensoft2/adm/etc/profile.d/modules.sh
module purge
module load samtools
module load bedtools
module load hisat2/2.0.1-beta
module load fastx_toolkit
module load Python/2.7.8
module load HTSeq/0.9.1

line=`head -n ${SAMPLE_NUM} ${SAMPLEDIR}/SampleFileFastq_tars.txt | tail -n 1`
SAMPLE=`echo $line | cut -d ' ' -f 1`
FASTQFILE_1=`echo $line | cut -d ' ' -f 2`
FASTQFILE_2=`echo $line | cut -d ' ' -f 3`
MERGE_REPLACE=`echo $line | cut -d ' ' -f 4`

TASKDIR="${TMPDIR}/${SAMPLE}"
mkdir ${TASKDIR}
OUTDIR="${HISATDIR}/Results/${SAMPLE}"
mkdir "${HISATDIR}/Results/${SAMPLE}"

FASTQFILE_NEW="${TASKDIR}/${SAMPLE}"

if [[ ${MERGE_REPLACE} == '' ]]; 
then 	
	echo 'copying'; 
	echo "cp ${FASTQFILE_1} ${FASTQFILE_NEW}.fq.gz"
	cp ${FASTQFILE_1} ${FASTQFILE_NEW}.fq.gz
elif [ ${MERGE_REPLACE} == '.replace' ];
then 
	echo 'replacing'; 
	echo "cp ${FASTQFILE_2} ${FASTQFILE_NEW}.fq.gz"
	cp ${FASTQFILE_2} ${FASTQFILE_NEW}.fq.gz
elif [ ${MERGE_REPLACE} == '.merge' ];
then 
	echo 'merging'; 
	echo "cp ${FASTQFILE_1} ${FASTQFILE_NEW}.fq.gz"
	echo "cat ${FASTQFILE_2} >> ${FASTQFILE_NEW}.fq.gz"
	cp ${FASTQFILE_1} ${FASTQFILE_NEW}.fq.gz
	cat ${FASTQFILE_2} >> ${FASTQFILE_NEW}.fq.gz
else echo 'problem'
fi

gunzip ${FASTQFILE_NEW}.fq.gz
# trim the last bp of each read (lower quality)

FASTQFILE_TRIMMED=${FASTQFILE_NEW}_trimmed.fq

fastx_trimmer -Q33 -t 1 -i ${FASTQFILE_NEW}.fq -o ${FASTQFILE_TRIMMED}
 rm ${FASTQFILE_NEW}.fq

ALIGNED="${TASKDIR}/hisat_aligned_${SAMPLE}"
UNALIGNED="${TASKDIR}/hisat_unmapped_${SAMPLE}.fq"
METRICS="${TASKDIR}/hisat_metrics_${SAMPLE}.txt"
KNOWNSPLICESITES="/pasteur/projets/policy01/evo_immuno_pop/Maxime/Splicing/HISAT2/test/splice_sites.txt"
INDEX="/pasteur/projets/policy01/evo_immuno_pop/Maxime/Splicing/HISAT2/grch37_snp/genome_snp"

######## ALIGN reads with Hisat2 (alignment done on the genome + snp)
hisat2 -x ${INDEX} -q -U ${FASTQFILE_TRIMMED} --known-splicesite-infile ${KNOWNSPLICESITES} --un ${UNALIGNED} --met-file ${METRICS} --dta | samtools view -Sb - > ${ALIGNED}.bam || exit 1
 rm ${FASTQFILE_TRIMMED}

samtools sort ${ALIGNED}.bam ${ALIGNED}.sort
samtools index ${ALIGNED}.sort.bam

GFF_EXONS='/pasteur/homes/mrotival/Annotation/All_Exons_merged_EnsV70.gff'
GFF_INTRONS='/pasteur/homes/mrotival/Annotation/All_Introns_merged_EnsV70.gff'

######## count Exonic Reads
htseq-count --stranded=no --format=bam ${ALIGNED}.sort.bam ${GFF_EXONS} --type=Exon --idattr=exon_id > ${OUTDIR}/${SAMPLE}_Exon_read_count.txt
######## count Intronic Reads
htseq-count --stranded=no --format=bam ${ALIGNED}.sort.bam ${GFF_INTRONS} --type=Intron --idattr=intron_id > ${OUTDIR}/${SAMPLE}_Intron_read_count.txt
######## count Gene Reads (aggregated)
htseq-count --stranded=no --format=bam ${ALIGNED}.sort.bam ${GFF_EXONS} --type=Exon --idattr=gene_id > ${OUTDIR}/${SAMPLE}_Gene_read_count.txt

 rm ${ALIGNED}.bam
 rm ${ALIGNED}.sort.bam
 rm ${ALIGNED}.sort.bam.bai
 rm ${UNALIGNED}

