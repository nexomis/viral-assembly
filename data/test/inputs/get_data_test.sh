#!/bin/bash

echo "get_data_test.sh"
echo "start : $(date)"

WD=data/test/inputs

rm -rf data/test/inputs/{fq,genome_ref,k2_HPRC_20230810,k2_standard_08gb_20240605,k2_viral_20240605,sars2_proteome.fa}

mkdir -p data/test/inputs/{fq,genome_ref,k2_HPRC_20230810,k2_standard_08gb_20240605,k2_viral_20240605}

echo "get fastq data"

# smpl1 : cleaned reads (host removed + trimmed?) : sars_cov2 - PE - fq.gz
srr="SRR10903401"
apptainer exec -W $WD/fq docker://ncbi/sra-tools:3.1.0 /bin/sh -c "\
  prefetch ${srr} && \
  fasterq-dump ${srr} && \
  rm -rf ${srr}/ && \
  gzip ${srr}*.fastq"

echo "compress using string"

apptainer exec -W $WD/fq docker://ghcr.io/nexomis/spring:1.1.1 spring -g --no-ids -q ill_bin -c -i ${srr}_1.fastq.gz ${srr}_1.fastq.gz -o s2.spring

# reference genome
echo "get reference genome"

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna.gz \
  -O $WD/genome_ref/sars_cov2_wuhan_refseq.fa.gz
gunzip $WD/genome_ref/sars_cov2_wuhan_refseq.fa.gz
cp $WD/genome_ref/sars_cov2_wuhan_refseq.fa $WD/genome_ref/vrl_genome_ref_s2.fa
gzip $WD/genome_ref/vrl_genome_ref_s2.fa

echo "get proteome"

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_protein.faa.gz \
  -O $WD/sars2_proteome.fa.gz
gunzip $WD/sars2_proteome.fa.gz
sed -i -r 's/YP_[0-9]+.[0-9] //g; s/ \[.*$//g ; s/ /_/g' $WD/sars2_proteome.fa

echo "get kraken2 databases"

for kurl in "https://zenodo.org/records/8339732/files/k2_HPRC_20230810.tar.gz" "https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20240605.tar.gz" "https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20240605.tar.gz"; do
  kfile="$(basename $kurl)"
  kdir="$(basename $kfile .tar.gz)"
  wget $kurl -O $WD/$kdir/$kfile
  tar -xvzf $WD/$kdir/$kfile -C $WD/$kdir/
  rm $WD/$kdir/$kfile
done
