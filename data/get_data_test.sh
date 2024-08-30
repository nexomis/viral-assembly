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
host_kraken_db_s1="k2_HPRC_20230810/db"
host_kraken_db_s2="k2_HPRC_20230810/db"

# dll kraken databases
mkdir k2_HPRC_20230810
cd k2_HPRC_20230810
wget https://zenodo.org/records/8339732/files/k2_HPRC_20230810.tar.gz
tar xvzf k2_HPRC_20230810.tar.gz
rm k2_HPRC_20230810.tar.gz
cd ..

mkdir k2_viral_20240605
cd k2_viral_20240605
wget https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20240605.tar.gz
tar xvzf k2_viral_20240605.tar.gz
rm k2_viral_20240605.tar.gz
cd ..

