

fq1=example_data/example_rna-seq_R1_001.fastq.gz
fq2=example_data/example_rna-seq_R2_001.fastq.gz
ncpu=2

mkdir bwa_ref
mkdir bwa_aln

cp example_out_0.2/context.fa bwa_ref/context.fa

# index
bwa index bwa_ref/context.fa

bwa aln -t ${ncpu} bwa_ref/context.fa ${fq1} > bwa_aln/`basename $fq1`.sai
bwa aln -t ${ncpu} bwa_ref/context.fa ${fq2} > bwa_aln/`basename $fq2`.sai
bwa sampe bwa_ref/context.fa bwa_aln/`basename $fq1`.sai bwa_aln/`basename $fq2`.sai ${fq1} ${fq2}  \
  | samtools view -uS - \
  | samtools sort  \
  > bwa_aln/example_rna-seq.bam
samtools index bwa_aln/example_rna-seq.bam


/code/STAR/2.6.1d/bin/Linux_x86_64_static/STAR \
  --outFileNamePrefix example_out_0.2/star/ \
  --limitOutSAMoneReadBytes 1000000 \
  --genomeDir example_out_0.2/STAR_idx \
  --readFilesCommand 'gzip -d -c -f' \
  --readFilesIn example_data/example_rna-seq_R1_001.fastq.gz example_data/example_rna-seq_R2_001.fastq.gz \
  --outSAMmode Full \
  --outFilterMultimapNmax -1 \
  --outSAMattributes Standard \
  --outSAMunmapped None \
  --outFilterMismatchNoverLmax 0.02 \
  -outFilterMismatchNoverReadLmax 0.02 \
  --outSAMtype BAM SortedByCoordinate \
  --runThreadN 12

