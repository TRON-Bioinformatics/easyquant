

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

