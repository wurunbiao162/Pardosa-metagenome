#!/usr/bin/env sh

#SBATCH -N 2 --ntasks 2 
#SBATCH --cpus-per-task=60
##SBATCH --nodelist=cluster02


mkdir -p result
mkdir -p seq 

#################### filter  ################
## The experimental design [metadata.txt] file is uploaded to the results folder
## The shotgun sequence data were uploaded to the [seq] folder.


## Quality control and assessment
time fastqc seq/*.fastq -t 20
multiqc -d seq/ -o result/qc


## Qaulity control and remove host contamination
db=~/db
mkdir -p ${db}/biobakery_workflows_databases/kneaddata_db_human_genome/
wget -c 
# build index
bowtie2-build Pardosa_pseudoannulata_genome.fa Pardosa_pseudoannulata_genome

for i in  `tail -n+2 result/metadata.txt | cut -f 1`;do \
kneaddata \
 --input rawdata/${i}.1.fastq.gz \
 --input rawdata/${i}.2.fastq.gz \
 --output seq \
 --threads 20 --output-prefix RemoveHost.${i} \
 --reference-db ${db}/biobakery_workflows_databases/kneaddata_db_human_genome \
 --serial --run-trf \
 --trimmomatic /data/home/wurunbiao/common/miniconda3/envs/metagenome/share/trimmomatic-0.39-1 \
 --remove-intermediate-output \
 --bowtie2-options "--very-sensitive --dovetail"
done

## Results will be in seq folder: such as RemoveHost.*.1.fastq.gz，RemoveHost.*.2.fastq.gz


# Rename the files
for i in  `tail -n+2 result/metadata.txt | cut -f 1`; do 
	mv seq/RemoveHost.${i}.1.fastq.gz seq/${i}_1.fastq
	mv seq/RemoveHost.${i}.2.fastq.gz seq/${i}_2.fastq
done

# Filtering out the human genome sequence
db=~/db
mkdir -p ${db}/BMTAGGER_INDEX
cd ${db}/BMTAGGER_INDEX
wget -c ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/*fa.gz
gunzip *fa.gz
cat *fa > hg38.fa
rm chr*.fa

bmtool -d hg38.fa -o hg38.bitmask
srprism mkindex -i hg38.fa -o hg38.srprism -M 100000
# Set the path to the database, "which config-metawrap" finds the location of the configuration file, and then you changed it such as 
# BMTAGGER_DB=/data/home/spider/db/BMTAGGER_INDEX

mkdir READ_QC
for F in seq/*_1.fastq; do R=${F%_*}_2.fastq; BASE=${F##*/}; SAMPLE=${BASE%_*}; metawrap read_qc -1 $F -2  $R -t 20 -o READ_QC/$SAMPLE & done
# Final sequence for decontamination: final_pure_reads_1.fastq，final_pure_reads_2.fastq.

## Move all final sequence files  to a new folder [CLEAN_READS]
mkdir CLEAN_READS
for i in READ_QC/*; do 
	b=${i#*/}
	mv ${i}/final_pure_reads_1.fastq CLEAN_READS/${b}_1.fastq
	mv ${i}/final_pure_reads_2.fastq CLEAN_READS/${b}_2.fastq
done



#################### Assembly-free metagenomic profiling  ################
conda install r kraken2 braken -y

conda activate kraken2
kraken2 --version # 2.1.2

## Kraken2 for classification of species
db=~/common/db
mkdir -p temp/kraken2

for i in `tail -n+2 result/metadata.txt | cut -f 1`; do \
kraken2 --db ${db}/MY_KRAKEN_DATABASE/minikraken_8GB_20200312 \
--paired temp/${i}_1.fastq temp/${i}_2.fastq \
--threads 40 --use-names --report-zero-counts \
--report temp/kraken2/${i}.report \
--output temp/kraken2/${i}.output
done


## krakentools changed report to mpa format.
conda install krakentools -c bioconda

for i in `tail -n+2 result/metadata.txt|cut -f1`
do
${db}/script/kreport2mpa.py -r temp/kraken2/${i}.report \
--display-header \
-o temp/kraken2/${i}.mpa 
done

mkdir -p result/kraken2
parallel -j 1 \
  'tail -n+2 temp/kraken2/{1}.mpa | sort | cut -f 2 | sed "1 s/^/{1}\n/" > temp/kraken2/{1}_count ' \
  ::: `tail -n+2 result/metadata.txt|cut -f1`

header=`tail -n 1 result/metadata.txt | cut -f 1`
echo $header
tail -n+2 temp/kraken2/${header}.mpa | sort | cut -f 1 | \
  sed "1 s/^/Taxonomy\n/" > temp/kraken2/0header_count
head -n3 temp/kraken2/0header_count

ls temp/kraken2/*count
paste temp/kraken2/*count > result/kraken2/tax_count.mpa
# Statistics
csvtk -t stat result/kraken2/tax_count.mpa


### Bracken for estimating abundance
db=~/common/db
tax=S
mkdir -p temp/bracken
for i in `tail -n+2 result/metadata.txt|cut -f 1`;do
bracken -d ${db}/MY_KRAKEN_DATABASE/minikraken_8GB_20200312 \
  -i temp/kraken2/${i}.report \
  -r 150 -l ${tax} -t 0 \
  -w temp/bracken/${i}.${tax}.bracken.report \
  -o temp/bracken/${i}
done


## The results are combined in a table.
parallel -j 1 \
  'tail -n+2 temp/bracken/{1}|sort|cut -f 6| sed "1 s/^/{1}\n/" > temp/bracken/{1}.count ' \
  ::: `tail -n+2 result/metadata.txt|cut -f1`

h=`tail -n1 result/metadata.txt|cut -f1`
tail -n+2 temp/bracken/${h}|sort|cut -f1 | \
  sed "1 s/^/Taxonomy\n/" > temp/bracken/0header.count

ls temp/bracken/*count | wc

paste temp/bracken/*count > result/kraken2/bracken.${tax}.txt

csvtk -t stat result/kraken2/bracken.${tax}.txt

Rscript ${db}/script/filter_feature_table.R \
  -i result/kraken2/bracken.${tax}.txt \
  -p 0.01 \
  -o result/kraken2/bracken.${tax}.0.01
# > 0.01(1%) species were only selected.
csvtk -t stat result/kraken2/bracken.${tax}.0.01

# Removal of chordates at phylum level
grep 'Chordata' result/kraken2/bracken.P.0.01
grep -v 'Chordata' result/kraken2/bracken.P.0.01 > result/kraken2/bracken.P.0.01-H

# Removal of P(Chordata) and S(Homo sapiens)
grep 'Homo sapiens' result/kraken2/bracken.${tax}.0.01
grep -v 'Homo sapiens' result/kraken2/bracken.${tax}.0.01 \
  > result/kraken2/bracken.${tax}.0.01-H

rm -rf temp/kraken2/*.output
# convert [bracken.S.norm] file to excel format [bracken.S.norm.xlsx].


## Taxonomy stackplot
# Phylum(P) and Species(S) level，results including 
# [bracken.P.stackplot.group.pdf],[bracken.P.stackplot.sample.pdf]
# [bracken.S.stackplot.group.pdf],[bracken.S.stackplot.sample.pdf]
sd=/data1/result/Pardosa_metag_project/script
tax=P
Rscript ${sd}/tax_stackplot.R \
  --input result/kraken2/bracken.${tax}.0.01-H --design result/metadata.txt \
  --group Group --output result/kraken2/bracken.${tax}.stackplot \
  --legend 10 --width 89 --height 79
  
tax=S
Rscript ${sd}/tax_stackplot.R \
  --input result/kraken2/bracken.${tax}.0.01-H --design result/metadata.txt \
  --group Group --output result/kraken2/bracken.${tax}.stackplot \
  --legend 10 --width 89 --height 79



## plot the diversity analysis and species composition 
cd /data1/result/Pardosa_metag_project
conda activate R 
sd=/data1/result/Pardosa_metag_project/script

### Alpha diversity
Rscript $sd/kraken2alpha.R -h

# --depth 0 means the minimum sample size
tax=S
Rscript $sd/otutab_rare.R \
  --input result/kraken2/bracken.${tax}.txt \
  --depth 0 --seed 1 \
  --normalize result/kraken2/bracken.${tax}.norm \
  --output result/kraken2/bracken.${tax}.alpha

# richness/chao1/ACE/shannon/simpson/invsimpson plot
for i in richness chao1 ACE shannon simpson invsimpson;do
Rscript $sd/alpha_boxplot.R -i result/kraken2/bracken.${tax}.alpha -a ${i} \
  -d result/metadata.txt 
  -n Group 
  -w 89 -e 59 \
  -o result/kraken2/${tax}
done

## Error in library(amplicon) : there is no package called ‘amplicon’
## To install amplicon package (https://github.com/microbiota/amplicon) and the server should be connected to the network.


### Beta diversity
mkdir -p kraken2/beta/
/c/db/win/usearch -beta_div result/kraken2/bracken.${tax}.norm \
-filename_prefix result/kraken2/beta/

## PCoA
# Available distances are bray_curtis, euclidean, jaccard, manhattan
dis=bray_curtis
Rscript $sd/beta_pcoa.R \
  --input result/kraken2/beta/${dis}.txt \
  --design result/metadata.txt \
  --group Group \
  --width 89 --height 59 \
  --output result/kraken2/pcoa.${dis}.pdf 



#################### Metagenome Assembly, Gene Prediction, and Annotation  ################
db=/data/home/spider/db
## MEGAHIT Assembly
    rm -rf temp/megahit
    cat CLEAN_READS/*_1.fastq > CLEAN_READS/left.fastq && cat CLEAN_READS/*_2.fastq > CLEAN_READS/right.fastq
    time megahit -t 60 \
    -1 CLEAN_READS/left.fastq \
    -2 CLEAN_READS/right.fastq  \ 
	-l 1000 --megahit \
    -o temp/megahit
	
    rm -rf  CLEAN_READS/left.fastq CLEAN_READS/right.fastq
    # the file name [final_assembly.fasta] was changed to [final.contigs.fa].
    mv temp/megahit/final_assembly.fasta temp/megahit/final.contigs.fa 
    ls -sh temp/megahit/final.contigs.fa

    seqkit stat temp/megahit/final.contigs.fa
    head -n6 temp/megahit/final.contigs.fa | cut -c1-60

    mkdir -p result/megahit/
    ln -f temp/megahit/final.contigs.fa result/megahit/
    # Delete temporary files
    rm -rf temp/megahit/intermediate_contigs/
    
### Gene prediction, cluster & quantitfy 
## metaProdigal for Gene prediction
    mkdir -p temp/prodigal

    time prodigal -i result/megahit/final.contigs.fa \
        -d temp/prodigal/gene.fa \
        -o temp/prodigal/gene.gff \
        -p meta -f gff > temp/prodigal/gene.log 2>&1 
    # check the log file
    tail temp/prodigal/gene.log
    # Counting the number of genes
    grep -c '>' temp/prodigal/gene.fa 
    grep -c 'partial=00' temp/prodigal/gene.fa 
    # Extraction of intact genes (complete fragments obtained with all genes intact, e.g. ring-forming bacterial genomes)
    grep 'partial=00' temp/prodigal/gene.fa | cut -f1 -d ' '| sed 's/>//' > temp/prodigal/full_length.id
    seqkit grep -f temp/prodigal/full_length.id temp/prodigal/gene.fa > temp/prodigal/full_length.fa
    seqkit stat temp/prodigal/full_length.fa


### Gene clustering/de-redundancy cd-hit  
    mkdir -p result/NR
    time cd-hit-est -i temp/prodigal/gene.fa \
        -o result/NR/nucleotide.fa \
        -aS 0.9 -c 0.95 -G 0 -g 0 -T 60 -M 0
		
    grep -c '>' result/NR/nucleotide.fa
    # Translated nucleic acids as corresponding protein sequences
    transeq -sequence result/NR/nucleotide.fa \
      -outseq result/NR/protein.fa -trim Y 

    sed -i 's/_1 / /' result/NR/protein.fa


### salmon   
    mkdir -p temp/salmon
    salmon -v # 1.4.0
    
    time salmon index \
      -t result/NR/nucleotide.fa \
      -p 60 \
      -i temp/salmon/index 
    
    time parallel -j 2 \
      "salmon quant \
        -i temp/salmon/index -l A -p 60 --meta \
        -1 CLEAN_READS/{1}_1.fastq \
        -2 CLEAN_READS/{1}_2.fastq \
        -o temp/salmon/{1}.quant" \
        ::: `tail -n+2 result/metadata.txt|cut -f1`
    
    # conbine
    mkdir -p result/salmon
    salmon quantmerge \
        --quants temp/salmon/*.quant \
        -o result/salmon/gene.TPM
    salmon quantmerge \
        --quants temp/salmon/*.quant \
        --column NumReads -o result/salmon/gene.count
    sed -i '1 s/.quant//g' result/salmon/gene.*
    
    # Preview results table
    head -n3 result/salmon/gene.*


### Functional gene annotation
## eggNOG(COG/KEGG/CAZy)
    # https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2
    
    mkdir -p temp/eggnog
    db=/data/home/spider/db
    time emapper.py --no_annot --no_file_comments --override \
      --data_dir ${db}/eggnog5 \
      -i result/NR/protein.fa \
      --cpu 60 -m diamond \
      -o temp/eggnog/protein

    time emapper.py \
      --annotate_hits_table temp/eggnog/protein.emapper.seed_orthologs \
      --data_dir ${db}/eggnog5 \
      --cpu 60 --no_file_comments --override \
      -o temp/eggnog/output

    sed '1 i Name\tortholog\tevalue\tscore\ttaxonomic\tprotein\tGO\tEC\tKO\tPathway\tModule\tReaction\trclass\tBRITE\tTC\tCAZy\tBiGG\ttax_scope\tOG\tbestOG\tCOG\tdescription' \
      temp/eggnog/output.emapper.annotations \
      > temp/eggnog/output

## Generating COG/KO/CAZy abundance summary tables
    mkdir -p result/eggnog
    python3 ${db}/script/summarizeAbundance.py -h
    python3 ${db}/script/summarizeAbundance.py \
      -i result/salmon/gene.TPM \
      -m temp/eggnog/output \
      -c '9,16,21' -s ',+,+*' -n raw \
      -o result/eggnog/eggnog
    sed -i 's/^ko://' result/eggnog/eggnog.KO.raw.txt
    # eggnog.CAZy.raw.txt  eggnog.COG.raw.txt  eggnog.KO.raw.txt
    

## GO analysis were performed with "GO analysis.r" in R software. 
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(stringr)

# go_files.txt file including query_name
egg <- read.table("C:/Users/Lenovo/Desktop/go_files.txt",sep="\t",header=T, stringsAsFactors = FALSE)
gene_ids <- egg$query_name
eggnog_lines_with_go <- egg$GOs!= ""

eggnog_annoations_go <- str_split(egg[eggnog_lines_with_go,]$GOs, ",")
gene_to_go <- data.frame(gene = rep(gene_ids[eggnog_lines_with_go], times = sapply(eggnog_annoations_go, length)), term = unlist(eggnog_annoations_go))
term2gene1 <- gene_to_go[, c(2, 1)]

term2gene <- buildGOmap(term2gene1)
go2term <- go2term(term2gene$GO)
go2ont <- go2ont(term2gene$GO)

# gene1.txt file including query_name and GOs.
gene1 <- read.table("C:/Users/Lenovo/Desktop/gene1.txt", header = FALSE, stringsAsFactors = FALSE)
gene1 <- gene1$V1[1:nrow(gene1)]
x <- enricher(gene = gene1, TERM2GENE = term2gene, TERM2NAME = go2term, pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.05)
#dotplot
dotplot(x, showCategory=20)


## KEGG analysis were performed on the website https://www.omicshare.com/
## Two documents need to be prepared in advance: beijingjiyinwenjian.txt file and mudijiyin.txt file
# 1.beijingjiyinwenjian.txt file including query and KEGG_ko.
# 2.mudijiyin.txt file including including query.




#################### Analysis of Genes Related to Carbohydrate and Antibiotic Resistance  ################
### Carbohydrate
    mkdir -p temp/dbcan2
    time diamond blastp \
      --db ${db}/dbCAN2/CAZyDB.07312020 \
      --query result/NR/protein.fa \
      --threads 60 --sensitive -e 1e-5 --outfmt 6 --max-target-seqs 1 --quiet \
    	--out temp/dbcan2/gene_diamond.f6
    wc -l temp/dbcan2/gene_diamond.f6

    mkdir -p result/dbcan2
    perl $db/script/format_dbcan2list.pl \
      -i temp/dbcan2/gene_diamond.f6 \
      -o temp/dbcan2/gene.list
    # Get abundance
    python3 ${db}/script/summarizeAbundance.py \
      -i result/salmon/gene.TPM \
      -m temp/dbcan2/gene.list \
      -c 2 -s ',' -n raw \
      -o result/dbcan2/TPM
    # spf format
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print a[$1],$0}' \
       ${db}/dbCAN2/CAZy_description.txt result/dbcan2/TPM.CAZy.raw.txt | \
      sed 's/^\t/Unannotated\t/' > result/dbcan2/TPM.CAZy.raw.spf


## Antibiotic Resistance
db=/data/home/wurunbiao/common/db
    mkdir -p temp/resfam
    time diamond blastp --db ${db}/protein/resfam/Resfams-proteins --query result/NR/protein.fa \
    	--outfmt 6 --threads 9 --max-target-seqs 1 --quiet \
    	--out temp/resfam/gene_diamond.f6
    
    mkdir -p result/resfam
    cut -f 1,2 temp/resfam/gene_diamond.f6 | uniq | \
      sed '1 i Name\tResGeneID' > temp/resfam/gene_fam.list

    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print a[$1],$0}' \
      temp/resfam/gene_fam.list result/salmon/gene.count | \
    	sed '/^\t/d' > result/resfam/resfam.count

    wc -l result/salmon/gene.count
    wc -l result/resfam/resfam.count # 369/652298=0.056%
    # spf format
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$4"\t"$3"\t"$2} NR>FNR{print a[$1],$0}' \
      ${db}/resfam/Resfams-proteins_class.tsv  result/resfam/resfam.count \
      > result/resfam/resfam.count.spf
# TPM.CAZy.raw.spf，resfam.count.spf and metadata.txt as input files to compare differences on STAMP software.



#################### Assembly metagenomic profiling  ################
# bin
	metawrap binning -o temp/INITIAL_BINNING -t 60 -a temp/megahit/final_assembly.fasta --metabat2 --maxbin2 --concoct CLEAN_READS/*.fastq

# Quality Control, Classification, and Annotation
	metawrap bin_refinement -o temp/BIN_REFINEMENT -t 60 \
		-A temp/INITIAL_BINNING/metabat2_bins/ \
		-B temp/INITIAL_BINNING/maxbin2_bins/ \
		-C temp/INITIAL_BINNING/concoct_bins/ \
		-c 50 -x 10 
	# Number of statistics
	cat temp/BIN_REFINEMENT/metawrap_50_10_bins.stats | awk '$2>50 && $3<10' | wc -l


# dRep for de-redundant species/strain genomes
    source ${soft}/bin/activate drep
	
    mkdir -p temp/DREP_IN
    ln -s `pwd`/temp/BIN_REFINEMENT/metawrap_50_10_bins/bin.* temp/DREP_IN/
    rename 'bin' 'mix_all' temp/DREP_IN/bin.*

    ls temp/DREP_IN/|cut -f 1 -d '_'|uniq -c
    ls temp/DREP_IN/|cut -f 2 -d '_'|cut -f 1 -d '.' |uniq -c

	# Redundancy by species level：
    mkdir -p temp/DREP95
    dRep dereplicate temp/DREP95/ \
      -g temp/DREP_IN/*.fa \
      -sa 0.95 -nc 0.30 -comp 50 -con 10 -p 3

# re-assembly bins
	source ~/wurunbiao/common/miniconda3/bin/activate
    source activate metawrap
	metawrap reassemble_bins -o temp/BIN_REASSEMBLY -1 CLEAN_READS/ALL_READS_1.fastq -2 CLEAN_READS/ALL_READS_2.fastq -t 60 -m 600 -c 70 -x 5 -b temp/DREP95/dereplicated_genomes/

# bin species annotation
	metawrap classify_bins -b temp/BIN_REASSEMBLY/reassembled_bins -o temp/BIN_CLASSIFICATION -t 60
	# BIN_CLASSIFICATION/bin_taxonomy.tab is the result file.

# bin function annotation
	metawrap annotate_bins -o temp/FUNCT_ANNOT -t 60 -b temp/BIN_REASSEMBLY/reassembled_bins/
##  FUNCT_ANNOT/bin_funct_annotations
##  FUNCT_ANNOT/bin_translated_genes,FUNCT_ANNOT/bin_untranslated_genes
##  PROKKA ouput FUNCT_ANNOT/prokka_out.


## assify  bins by GTDB-tk
    source ~/wurunbiao/common/miniconda3/bin/activate 
    source activate gtdbtk
	
    mkdir -p temp/gtdb_classify
    # 10个基因组，24p，100min 152 G内存
    gtdbtk classify_wf \
        --genome_dir temp/BIN_REASSEMBLY/reassembled_bins \
        --out_dir temp/gtdb_classify \
        --extension fa \
        --prefix tax \
        --cpus 24



#################### Phylogenetic Analysis  ################
    mkdir -p temp/gtdb_infer
    gtdbtk infer \
        --msa_file temp/gtdb_classify/tax.bac120.user_msa.fasta \
        --out_dir temp/gtdb_infer \
        --prefix tax \
        --cpus 24

## table2itol
# tax.bac120.summary.tsv, Widb.csv
    mkdir -p result/itol
    tail -n+2 temp/gtdb_classify/tax.bac120.summary.tsv|cut -f 1-2|sed 's/;/\t/g'|sed '1 s/^/ID\tDomain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n/' \
      > result/itol/tax.txt

    sed 's/,/\t/g;s/.fa//' temp/DREP95/data_tables/Widb.csv|cut -f 1-7,11|sed '1 s/genome/ID/' \
      > result/itol/genome.txt

    awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print $0,a[$1]}' result/itol/genome.txt result/itol/tax.txt|cut -f 1-8,10- > result/itol/annotation.txt

    cd result/itol/
    db=/data/home/spider/db
    Rscript ${db}/script/table2itol.R -a -c double -D plan1 -i ID -l Species -t %s -w 0.5 annotation.txt
