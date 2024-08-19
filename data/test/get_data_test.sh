cd data/test/
mkdir fq/ genome_ref/

# smpl1 : cleaned reads (host removed + trimmed?) : sars_cov2 - PE - fq.gz
srr="SRR10903401"
docker run -u $(id -u):$(id -g) -v $PWD/fq/:$PWD/fq/ -w $PWD/fq/ ncbi/sra-tools:3.1.0 /bin/sh -c "\
  prefetch ${srr} && \
  fasterq-dump ${srr} && \
  rm -rf ${srr}/ && \
  gzip ${srr}*.fastq"

# smpl2 : copie of smpl1
cp fq/${srr}_1.fastq.gz fq/s2_R1.fq.gz 
cp fq/${srr}_2.fastq.gz fq/s2_R2.fq.gz 

# reference genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna.gz \
  -O genome_ref/sars_cov2_wuhan_refseq.fa.gz
cp genome_ref/vrl_genome_ref_s2.fa.gz
gunzip genome_ref/sars_cov2_wuhan_refseq.fa.gz genome_ref/vrl_genome_ref_s2.fa.gz

# samplesheet.csv
echo -e "\
sample_name,path_fq_R1,path_fq_R2,path_ref_genome,taxid_host_genome,assembly_type\n\
s1,data/test/fq/SRR10903401_1.fastq.gz,data/test/fq/SRR10903401_2.fastq.gz,data/test/genome_ref/sars_cov2_wuhan_refseq.fa,9606,corona\n\
s2,data/test/fq/s2_R1.fq.gz,data/test/fq/s2_R2.fq.gz,data/test/genome_ref/vrl_genome_ref_s2.fa,9606,rnaviral\
" > samplesheet.csv

# exec:
 git/viral-assembly$ NXF_VER=24.04.4 nextflow run main.nf \
  --input data/test/samplesheet.csv \
  --skip_primary \
  --mem_high 11.GB \
  -profile docker
