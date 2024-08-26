mkdir test/
mkdir test/fq/ test/genome_ref/

# smpl1 : cleaned reads (host removed + trimmed?) : sars_cov2 - PE - fq.gz
srr="SRR10903401"
docker run -u $(id -u):$(id -g) -v $PWD/test/fq/:$PWD/test/fq/ -w $PWD/test/fq/ ncbi/sra-tools:3.1.0 /bin/sh -c "\
  prefetch ${srr} && \
  fasterq-dump ${srr} && \
  rm -rf ${srr}/ && \
  gzip ${srr}*.fastq"

# smpl2 : cp of smpl1
cp test/fq/${srr}_1.fastq.gz test/fq/s2_R1.fq.gz
cp test/fq/${srr}_2.fastq.gz test/fq/s2_R2.fq.gz

# reference genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna.gz \
  -O test/genome_ref/sars_cov2_wuhan_refseq.fa.gz
gunzip test/genome_ref/sars_cov2_wuhan_refseq.fa.gz
cp test/genome_ref/sars_cov2_wuhan_refseq.fa test/genome_ref/vrl_genome_ref_s2.fa

# samplesheet.csv
host_kraken_db_s1="/path/to/host_kraken_db/s1/"
host_kraken_db_s2="/path/to/host_kraken_db/s2/"
echo -e "\
sample_name,path_fq_R1,path_fq_R2,path_ref_genome,host_kraken_db,assembly_type\n\
s1,test/fq/SRR10903401_1.fastq.gz,test/fq/SRR10903401_2.fastq.gz,test/genome_ref/sars_cov2_wuhan_refseq.fa,${host_kraken_db_s1},corona\n\
s2,test/fq/s2_R1.fq.gz,test/fq/s2_R2.fq.gz,test/genome_ref/vrl_genome_ref_s2.fa,${host_kraken_db_s2},rnaviral\
" > test/samplesheet.csv
