
star_bin="/code/STAR/2.6.1d/bin/Linux_x86_64_static/STAR"
samtools_bin="/code/samtools/1.9/samtools"

fq1=example_data/example_rna-seq_R1_001.fastq.gz
fq2=example_data/example_rna-seq_R2_001.fastq.gz
ncpu=12
aln_out_dir="star_aln"

mkdir star_ref
mkdir ${aln_out_dir}

cp example_out_0.2/context.fa star_ref/context.fa

# Generate index
${star_bin} --runMode genomeGenerate \
  --runThreadN ${ncpu} \
  --genomeSAindexNbases 4 \
  --genomeDir star_ref/ \
  --genomeFastaFiles star_ref/context.fa

# map reads
# for STAR manual see: https://raw.githubusercontent.com/alexdobin/STAR/master/doc/STARmanual.pdf
${star_bin} \
  --outFileNamePrefix ${aln_out_dir}/ \
  --limitOutSAMoneReadBytes 1000000 \
  --genomeDir star_ref \
  --readFilesCommand 'gzip -d -c -f' \
  --readFilesIn ${fq1} ${fq2} \
  --outSAMmode Full \
  --outFilterMultimapNmax -1 \
  --outSAMattributes NH HI AS nM NM \
  --outSAMunmapped None \
  --outFilterMismatchNoverLmax 0.02 \
  --outFilterMismatchNoverReadLmax 0.02 \
  --outSAMtype BAM SortedByCoordinate \
  --runThreadN ${ncpu}

# index output bam
${samtools_bin} index ${aln_out_dir}/Aligned.sortedByCoord.out.bam

