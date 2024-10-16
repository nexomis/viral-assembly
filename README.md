# viral-assembly

## Description

Assemble viral genome and perform the associated QC. Starting from raw reads and one pre-built reference database in standard format (cf. README.md).

![Schema not visible](./schema.drawio.svg)

## Test with simple data set
First, edit manually 'data/get_data_test.sh': replace path 'host_kraken_db_s1' and 'host_kraken_db_s2'.

````
cd data/
bash get_data_test.sh

kraken_db_qc="/path/to/kraken/db/for/QC/"
NXF_VER=23.10.1 nextflow run ../main.nf \
  -profile docker \
  --input test/samplesheet.csv \
  --mem_high 8.GB \
  --kraken2_db k2_viral_20240605 \
  --output_dir test/out_dir/
```

## Citations

TODO add citations

