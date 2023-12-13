# fasta file
wget -O data/output/spliceai/hg38.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip data/output/spliceai/hg38.fa.gz

#chromHMM 15 core state
tar -zxvf data/annotation_database/master38.chromhmm.bedg.tar.gz

# AnnoVar annotation file
wget -O data/output/annovar/humandb/hg38_refGene.txt.gz http://www.openbioinformatics.org/annovar/download/hg38_refGene.txt.gz
wget -O data/output/annovar/humandb/hg38_refGeneMrna.fa.gz http://www.openbioinformatics.org/annovar/download/hg38_refGeneMrna.fa.gz
wget -O data/output/annovar/humandb/hg38_refGeneVersion.txt.gz http://www.openbioinformatics.org/annovar/download/hg38_refGeneVersion.txt.gz
wget -O data/output/annovar/humandb/hg38_dbnsfp33a.txt.gz http://www.openbioinformatics.org/annovar/download/hg38_dbnsfp33a.txt.gz
wget -O data/output/annovar/humandb/hg38_dbnsfp33a.txt.idx.gz http://www.openbioinformatics.org/annovar/download/hg38_dbnsfp33a.txt.idx.gz
wget -O data/output/annovar/humandb/hg38_dbnsfp41a.txt.gz http://www.openbioinformatics.org/annovar/download/hg38_dbnsfp41a.txt.gz
wget -O data/output/annovar/humandb/hg38_dbnsfp41a.txt.idx.gz http://www.openbioinformatics.org/annovar/download/hg38_dbnsfp41a.txt.idx.gz
wget -O data/output/annovar/humandb/hg38_gnomad30_genome.txt.gz http://www.openbioinformatics.org/annovar/download/hg38_gnomad30_genome.txt.gz
wget -O data/output/annovar/humandb/hg38_gnomad30_genome.txt.idx.gz http://www.openbioinformatics.org/annovar/download/hg38_gnomad30_genome.txt.idx.gz
gunzip data/output/annovar/humandb/hg38_refGene.txt.gz
gunzip data/output/annovar/humandb/hg38_refGeneMrna.fa.gz
gunzip data/output/annovar/humandb/hg38_refGeneVersion.txt.gz
gunzip data/output/annovar/humandb/hg38_dbnsfp33a.txt.gz
gunzip data/output/annovar/humandb/hg38_dbnsfp33a.txt.idx.gz
gunzip data/output/annovar/humandb/hg38_dbnsfp41a.txt.gz
gunzip data/output/annovar/humandb/hg38_dbnsfp41a.txt.idx.gz
gunzip data/output/annovar/humandb/hg38_gnomad30_genome.txt.gz
gunzip data/output/annovar/humandb/hg38_gnomad30_genome.txt.idx.gz