## Config
```shell
#ref cel-seq pipeline
conda install biopython
conda install -c bioconda bowtie2
conda install pysam
```
## Batch fetch raw data
```shell
#BioProject PRJNA272543
cat filereport_read_run_PRJNA272543_tsv.txt |awk 'NR>1{print $(NF-4)}'>sra.url

cat filereport_read_run_PRJNA272543_tsv.txt |awk 'NR>1{print $(NF-5)}'>fqftp.url

outputdir=./database/project/2016-HE-EmbryoDevloment/cel-seq/126-seq/rawdata/rawdata-126

cat  fq.url |while read id 
do 
echo "ascp -k 1 -QT -l 300m -p3301  \
-i /home/zk/anaconda3/etc/asperaweb_id_dsa.openssh \
era-fasp@${id} ${outputdir}"
done >fq.download.sh

#
cat fq.txt |while read id 
do
ascp -QT -l 300m -P33001  \
-i ~/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh   \
era-fasp@$id  .
done

#md5
awk 'NR>1{print$5"\t"$3}' filereport_read_run_PRJNA272543_tsv.txt > md5.txt
less -S md5.txt
md5sum -c md5.txt

```
## Fastqc
```shell
#
qcdir=./database/project/2016-HE-EmbryoDevloment/cel-seq/126-seq/rawdata/fastqc
fqdir=./database/project/2016-HE-EmbryoDevloment/cel-seq/126-seq/rawdata/
fastqc -t 2 -o $qcdir $fqdir/SRR*.fastq.gz
multiqc -o ./ SRR*.zip
```
## Bowtie2
```shell
cd ./database/genome/NCBI/HE/2017
bowtie2-build nHd_3.1_cds_from_spikein_genomic.fna bcs_Hd
ls *.gz | cut -c 8-10 | sort -u | while read id; do \
bowtie2  -x ./database/genome/NCBI/HE/2017/bowtiecs/bcs_Hd \
-U SRR1755${id}.fastq.gz \
-S ./database/project/2016-HE-EmbryoDevloment/126-seq/bowtiecs/SRR1755${id}.sam \
-p 3  2>SRR1755${id}Align.summary &done
```
