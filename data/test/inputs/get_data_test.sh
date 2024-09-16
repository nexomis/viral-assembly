mkdir -p data/test/inputs

cd data/test/inputs

mkdir -p fq

# smpl1 : cleaned reads (host removed + trimmed?) : sars_cov2 - PE - fq.gz
srr="SRR10903401"
docker run -u $(id -u):$(id -g) -v $PWD/fq/:$PWD/fq/ -w $PWD/fq/ ncbi/sra-tools:3.1.0 /bin/sh -c "\
  prefetch ${srr} && \
  fasterq-dump ${srr} && \
  rm -rf ${srr}/ && \
  gzip ${srr}*.fastq"

docker run -it -u $UID:$GID -w $PWD/fq -v $PWD/fq:$PWD/fq ghcr.io/nexomis/spring:1.1.1 spring -g --no-ids -q ill_bin -c -i ${srr}_1.fastq.gz ${srr}_1.fastq.gz -o s2.spring

# reference genome
mkdir -p genome_ref

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna.gz \
  -O genome_ref/sars_cov2_wuhan_refseq.fa.gz
gunzip genome_ref/sars_cov2_wuhan_refseq.fa.gz
cp genome_ref/sars_cov2_wuhan_refseq.fa genome_ref/vrl_genome_ref_s2.fa

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

mkdir k2_standard_08gb_20240605
cd k2_standard_08gb_20240605
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20240605.tar.gz
tar xvzf k2_standard_08gb_20240605.tar.gz
rm k2_standard_08gb_20240605.tar.gz
cd ..

cd ../../..
