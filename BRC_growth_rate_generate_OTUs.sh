#!/bin/bash
#
#SBATCH --job-name Nimrod_shivta_site_analysis
#SBATCH --cpus-per-task=1
#SBATCH --mem=4000
#SBATCH --output=Shivta-%j.out
#SBATCH --error=Shivta-%j.err

#sbatch Nimrod_shivta_site_analysis.sh $DATASET $PLATFORM $DATAFILETYPE

#. /etc/profile

#title          :Nimrod_shivta_site_analysis.sh
#description    :Analysis of Nimrod_shivta_site samples
#author         :Roey Angel
#date           :20160315
#version        :1.0    
#usage          :sbatch Nimrod_shivta_site_analysis.sh ${DATASET} ${PLATFORM} ${DATAFILETYPE}
#notes          :       
#bash_version   :4.3.30(1)-release
#============================================================================

module load usearch
module load mothur
module load qiime

. /etc/profile

# Set parameters
#DATASET=$1
DATASET="Shivta_site"
#PLATFORM=$2
PLATFORM="MiSeq"
#DATAFILETYPE=$3
DATAFILETYPE="FASTAQ"
# MINAMPLENGTH=$4
MINAMPLENGTH="400"

declare -a SAMPLES
SAMPLES=("SlpLime1S1" "SlpLime2S2" "SlpClk1S3" "SlpClk2S4" "SlpClk3S5" "SlpClk4S6" "SlpClk5S7" "SlpClk6S8" "SlpClk7S9" "SlpClk8S10" "SlpClk9S11" "ShivSClk1S12" "ShivSLimek1S13" "ShivSLimek2S14" "ShivSClk2S15" "ShivSClk3S16" "ShivSClk4S17" "ShivSClk5S18" "ShivSLimek3S19" "ShivSLimek4S20" "ShivSLimek5S21" "SlpLime3S22" "SlpLime4S23" "SlpLime5S24" "ShivSLimek6S25")
COUNTMIN="2"
OTUdef="97"
OTU_id=`echo -e "(100-${OTUdef})/100" | bc`
OTU_radius=`echo -e "100-${OTUdef}" | bc`
CLASSIFY="Y"
CHIMCHECK="Y"

# Set paths
MOTHUR=/apps/mothur/1.36.1/
USEARCH=/apps/usearch/8.1.1861/
DBS=~/DBs/		       	
SCRIPTS=~/Scripts/MiSeq/
# RUNDIR=/scratch/angel/Hyperarid_erosion/shivta_site_samples/
RUNDIR=`pwd`"/"
WORKDIR=/tmp/${DATASET}_work/
RESULTSDIR=/tmp/${DATASET}_results/

# Create dirs
if [[ -d ${WORKDIR} ]]
then
 rm -r ${WORKDIR}
fi
if [[ -d ${RESULTSDIR} ]]
then
 rm -r ${RESULTSDIR}
fi

mkdir ${WORKDIR}
mkdir ${RESULTSDIR}

# Process raw files according to platform and data type
echo -e "Analysis log of the ${DATASET} dataset.\nPlatform=${PLATFORM}, filetype=${DATAFILETYPE}" > ${RUNDIR}${DATASET}.log
echo -e "===========================================\n" >> ${RUNDIR}${DATASET}.log
echo -e "Using: " >> ${RUNDIR}${DATASET}.log
${USEARCH}usearch >> ${RUNDIR}${DATASET}.log
echo -e "Cite: Edgar, 2013; Edgar & Flyvbjerg, 2015" >> ${RUNDIR}${DATASET}.log
echo -e "-- \n" >> ${RUNDIR}${DATASET}.log

for fastqfile in *.fastq
do
    gzip ${fastqfile}
done

cp *.gz ${WORKDIR}
cd ${WORKDIR}
gunzip ${WORKDIR}*.gz 

echo -e "The following files were found in the folder and will be analysed:" >> ${RUNDIR}${DATASET}.log
FASTQFILES=`ls *.fastq`
echo -e $FASTQFILES >> ${RUNDIR}${DATASET}.log
echo -e "-- \n" >> ${RUNDIR}${DATASET}.log
  
# Merge fastq files according to run and quality filter
echo -e "Merge fastq pairs by samples into one file" >> ${RUNDIR}${DATASET}.log
echo -e "------------------------------------------" >> ${RUNDIR}${DATASET}.log

echo -e "Prepending sample names to file names:" >> ${RUNDIR}${DATASET}.log
FASTQFILESR1=($(ls *R1*.fastq)) # this is bad programming but what the hell
FASTQFILESR2=($(ls *R2*.fastq))
for (( i=0; i<${#SAMPLES[@]}; i++ )); do
    mv ${FASTQFILESR1[i]} ${SAMPLES[i]}${FASTQFILESR1[i]}
    mv ${FASTQFILESR2[i]} ${SAMPLES[i]}${FASTQFILESR2[i]}
done
echo -e "The following samples names have been prepended: ${SAMPLES[@]}\n" >> ${RUNDIR}${DATASET}.log

echo -e "Merge fastq pairs:" >> ${RUNDIR}${DATASET}.log
#echo -e "usearch -fastq_mergepairs *_R1_*.fastq -relabel @  -report mergereport.txt -fastqout ${DATASET}_merged.fq" >> ${RUNDIR}${DATASET}.log
#${USEARCH}usearch -fastq_mergepairs *_R1_*.fastq -relabel @ -report mergereport.txt -fastqout ${DATASET}_merged.fq
#cat mergereport.txt >> ${RUNDIR}${DATASET}.log
module load qiime

#multiple_split_libraries_fastq.py -i ${WORKDIR} -o trimmed -p ${RUNDIR}param.file.txt
echo -e "multiple_join_paired_ends.py -i ${WORKDIR} -o joined -p ${RUNDIR}param.file.txt" >> ${RUNDIR}${DATASET}.log
multiple_join_paired_ends.py -i ${WORKDIR} -o joined -p ${RUNDIR}param.file.txt 
module unload qiime

mkdir ./joined/LowQualMerged
for dir in ./joined/*/
do
    dir=${dir%*/}
    cd ${dir}
    currentDir=`echo ${dir##*/} | awk -F'[_.]' '{print $1}'`
    mv fastqjoin.join.fastq ${currentDir}.fastq
    echo -e "\nCharacters distribuion in ${currentDir}.fastq:" >> ${RUNDIR}${DATASET}.log
    echo -e "usearch -fastq_chars  ${currentDir}.fastq -log chars.log"
    ${USEARCH}usearch -fastq_chars  ${currentDir}.fastq -log chars.log
    cat chars.log >> ${RUNDIR}${DATASET}.log
    
    echo -e "\nQuality statistics:" >> ${RUNDIR}${DATASET}.log
    echo -e "usearch -fastq_eestats2 ${DATASET}_merged.fq -output eestats2.txt -length_cutoffs 200,300,10" >> ${RUNDIR}${DATASET}.log
    ${USEARCH}usearch -fastq_eestats2 ${DATASET}_merged.fq -output eestats2.txt -length_cutoffs 200,300,10
    cat eestats2.txt >> ${RUNDIR}${DATASET}.log
    ${USEARCH}usearch -fastq_filter ${currentDir}.fastq  -relabel @ -fastq_minlen ${MINAMPLENGTH} -fastq_maxee 10 -fastaout ../${currentDir}.fa
    ${USEARCH}usearch -fastq_filter ${currentDir}.fastq  -relabel @ -fastq_minlen ${MINAMPLEGTH}  -fastaout ../LowQualMerged/${currentDir}.fa
    
    cd ${WORKDIR}
done

echo -e "-- \n" >> ${RUNDIR}${DATASET}.log
cd ${WORKDIR}
cat ./joined/*.fa >> ${DATASET}_filtered.fa
cat ./joined/LowQualMerged/*.fa >> ${DATASET}_LQ_filtered.fa

echo -e "Quality statistics" >> ${RUNDIR}${DATASET}.log
${USEARCH}usearch -fastq_eestats2 ${DATASET}_filtered.fa -output eestats2.txt -length_cutoffs 200,300,10
echo -e "usearch -fastq_eestats2 ${DATASET}_filtered.fa -output eestats2.txt -length_cutoffs 200,300,10" >> ${RUNDIR}${DATASET}.log
cat eestats2.txt >> ${RUNDIR}${DATASET}.log

${MOTHUR}mothur "#summary.seqs(fasta=${DATASET}_filtered.fa, processors=10)" > tmp.log
grep -m 1 'mothur >' tmp.log >> ${RUNDIR}${DATASET}.log 
grep -m 1 -A1000 'Start' tmp.log | head -n -3  >> ${RUNDIR}${DATASET}.log
echo -e "-- \n" >> ${RUNDIR}${DATASET}.log

# Dereplicate sequences for OTU table
echo -e "Dereplicate sequences for OTU table" >> ${RUNDIR}${DATASET}.log
echo -e "-----------------------------------" >> ${RUNDIR}${DATASET}.log
echo -e "usearch -derep_fulllength ${DATASET}_filtered.fa -minuniquesize ${COUNTMIN} -fastaout ${DATASET}_unique.fa -sizeout --log tmp.log" >> ${RUNDIR}${DATASET}.log
${USEARCH}usearch -derep_fulllength ${DATASET}_filtered.fa -minuniquesize ${COUNTMIN} -fastaout ${DATASET}_unique.fa -sizeout --log tmp.log
if [[ ${COUNTMIN} -gt 1 ]]; then
  echo -e "Sequences smaller than ${COUNTMIN} will not be allowed to form OTUs" >> ${RUNDIR}${DATASET}.log
fi
cat tmp.log >> ${RUNDIR}${DATASET}.log
echo -e "-- \n" >> ${RUNDIR}${DATASET}.log

# Sort sequences by size
echo -e "Sort sequences by size" >> ${RUNDIR}${DATASET}.log
echo -e "----------------------" >> ${RUNDIR}${DATASET}.log
echo -e "usearch -sortbysize ${DATASET}_unique.fa -fastaout ${DATASET}_seqs_sorted.fa --log tmp.log" >> ${RUNDIR}${DATASET}.log
${USEARCH}usearch -sortbysize ${DATASET}_unique.fa -fastaout ${DATASET}_seqs_sorted.fa --log tmp.log
cat tmp.log >> ${RUNDIR}${DATASET}.log
echo -e "-- \n" >> ${RUNDIR}${DATASET}.log

# Find OTU representatives
echo -e "Find OTU representatives" >> ${RUNDIR}${DATASET}.log
# if chimeras are being filtered:
if [[ $CHIMCHECK == "Y" ]]
then
 echo -e "-----------------------" >> ${RUNDIR}${DATASET}.log
 echo -e "Chimeras are being checked" >> ${RUNDIR}${DATASET}.log
 echo -e " usearch -cluster_otus ${DATASET}_seqs_sorted.fa -otu_radius_pct $OTU_radius -sizeout -otus ${DATASET}_otus.fa -relabel OTU -maxrejects 0 -maxaccepts 0 --log tmp.log" >> ${RUNDIR}${DATASET}.log
 ${USEARCH}usearch -cluster_otus ${DATASET}_seqs_sorted.fa -otu_radius_pct $OTU_radius -sizeout -otus ${DATASET}_otus.fa -relabel OTU -maxrejects 0 -maxaccepts 0 --log tmp.log
 cat tmp.log >> ${RUNDIR}${DATASET}.log
 echo -e "\n" >> ${RUNDIR}${DATASET}.log
 # make a ref-based chimera check
 echo -e "Make a ref-based chimera check:" >> ${RUNDIR}${DATASET}.log
 echo e- "usearch -uchime_ref ${DATASET}_otus.fa -db ${DBS}gold.fa -nonchimeras nochimeras.fa -strand plus --log tmp.log" >> ${RUNDIR}${DATASET}.log
 ${USEARCH}usearch -uchime_ref ${DATASET}_otus.fa -db ${DBS}gold.fa -nonchimeras ${DATASET}_otuReps.fa -strand plus --log tmp.log
 cat tmp.log >> ${RUNDIR}${DATASET}.log
 echo -e "-- \n" >> ${RUNDIR}${DATASET}.log
 
else
 echo -e "Chimeras are not being checked" >> ${RUNDIR}${DATASET}.log
 echo -e "------------------------------" >> ${RUNDIR}${DATASET}.log
 echo -e "usearch -cluster_smallmem ${DATASET}_seqs_sorted.fa -id $OTU_id -sizeout -gapopen 10.0I/10.0E -gapext 1.0I/1.0E -centroids ${DATASET}_otus.fa -usersort -maxrejects 0 -maxaccepts 0 --log tmp.log" >> ${RUNDIR}${DATASET}.log
 ${USEARCH}usearch -cluster_smallmem ${DATASET}_seqs_sorted.fa -id $OTU_id -sizeout -gapopen 10.0I/10.0E -gapext 1.0I/1.0E -centroids ${DATASET}_otuReps.fa -usersort -maxrejects 0 -maxaccepts 0 --log tmp.log
 cat tmp.log >> ${RUNDIR}${DATASET}.log
 echo -e "-- \n" >> ${RUNDIR}${DATASET}.log
fi
cp ${DATASET}_otuReps.fa ${RESULTSDIR}${DATASET}_otuReps.fa

# Determine abundance of OTUs in samples  
echo -e "Determine abundance of OTUs in samples" >> ${RUNDIR}${DATASET}.log
echo -e "--------------------------------------" >> ${RUNDIR}${DATASET}.log
echo -e "usearch -usearch_global ${DATASET}_filtered.fa -db ${RESULTSDIR}${DATASET}_otuReps.fa -strand plus -id $OTU_id -maxrejects 0 -maxaccepts 0 -maxhits 1 -uc ${DATASET}_readmap.uc -matched ${RESULTSDIR}${DATASET}_matched_seqs.fa -mothur_shared_out ${DATASET}_otuTab.shared -otutabout ${DATASET}_otuTab.txt" >> ${RUNDIR}${DATASET}.log
${USEARCH}usearch -usearch_global ${DATASET}_filtered.fa -db ${RESULTSDIR}${DATASET}_otuReps.fa -strand plus -id $OTU_id -maxrejects 0 -maxaccepts 0 -maxhits 1 -uc ${DATASET}_readmap.uc -matched ${RESULTSDIR}${DATASET}_matched_seqs.fa -mothur_shared_out ${DATASET}_otuTab.shared -otutabout ${DATASET}_otuTab.txt
cp ${DATASET}_otuTab.* ${RESULTSDIR}
echo -e `grep -v -c "*" ${DATASET}_readmap.uc`" merged read pairs map to clusters" >> ${RUNDIR}${DATASET}.log
echo -e "-- \n" >> ${RUNDIR}${DATASET}.log

# How many reads did not match any OTU?
echo -e "How many sequences did not match any OTU?" >> ${RUNDIR}${DATASET}.log
echo -e "-----------------------------------------" >> ${RUNDIR}${DATASET}.log
grep -c "^N" ${DATASET}_readmap.uc >> ${RUNDIR}${DATASET}.log
echo -e "-- \n" >> ${RUNDIR}${DATASET}.log

# Classify sequences 
echo -e "Classify sequences" >> ${RUNDIR}${DATASET}.log
echo -e "------------------" >> ${RUNDIR}${DATASET}.log
if [[ $CLASSIFY == "Y" ]]
then
 echo -e "Classifying OTUs:" >> ${RUNDIR}${DATASET}.log
 cp ${DBS}silva.nr_v119_classifier.zip ./
 unzip silva.nr_v119_classifier.zip
 cp ${RESULTSDIR}/${DATASET}_otuReps.fa ${DATASET}_otuReps.fa
 ${MOTHUR}mothur "#classify.seqs(fasta=${DATASET}_otuReps.fa,template=silva.nr_v119.align, taxonomy=silva.nr_v119.tax, processors=10)" > tmp.log
 grep -m 1 'mothur >' tmp.log >> ${RUNDIR}${DATASET}.log 
 grep -m 1 -A1000 'Output File Names:' tmp.log | head -n -3  >> ${RUNDIR}${DATASET}.log   
 echo -e "-- \n" >> ${RUNDIR}${DATASET}.log
 cp ${DATASET}_otuReps.nr_v119.wang.taxonomy ${RESULTSDIR}${DATASET}_silva.nrv119.taxonomy
fi

# Creat singleton table
echo -e "Creating an OTU table of rejected singletons" >> ${RUNDIR}${DATASET}.log
echo -e "----------------------------------------------" >> ${RUNDIR}${DATASET}.log
echo -e "This is used for evaluation only. It should worry us if the singletons show some pattern." >> ${RUNDIR}${DATASET}.log
# Dereplicate sequences for singleton table
echo -e "Dereplicate sequences for singleton table:" >> ${RUNDIR}${DATASET}.log
echo -e "{USEARCH}usearch -derep_fulllength ${RESULTSDIR}${DATASET}_matched_seqs.fa -minuniquesize 1 -fastaout ${DATASET}_unique.fa -sizeout --log tmp.log" >> ${RUNDIR}${DATASET}.log
${USEARCH}usearch -derep_fulllength ${RESULTSDIR}${DATASET}_matched_seqs.fa -minuniquesize 1 -fastaout ${DATASET}_unique.fa -sizeout --log tmp.log

cat tmp.log >> ${RUNDIR}${DATASET}.log
echo -e "-- \n" >> ${RUNDIR}${DATASET}.log

cp ${DATASET}_unique.fa ${RESULTSDIR}${DATASET}_unique.fa

# Determine abundance of singletons in samples  
echo -e "Determine abundance of singletons in samples:" >> ${RUNDIR}${DATASET}.log
echo -e "--------------------------------------------" >> ${RUNDIR}${DATASET}.log
echo -e "usearch -usearch_global ${RESULTSDIR}${DATASET}_matched_seqs.fa -db ${RESULTSDIR}${DATASET}_unique.fa -strand plus -id 1.0 -gapopen 10.0I/10.0E -gapext 1.0I/1.0E -maxrejects 0 -maxhits 1 -uc ${DATASET}_unique.uc -mothur_shared_out ${DATASET}_uniqueOtuTab.shared -otutabout ${DATASET}_uniqueOtuTab.txt" >> ${RUNDIR}${DATASET}.log
${USEARCH}usearch -usearch_global ${RESULTSDIR}${DATASET}_matched_seqs.fa -db ${RESULTSDIR}${DATASET}_unique.fa -strand plus -id 1.0 -gapopen 10.0I/10.0E -gapext 1.0I/1.0E -maxrejects 0 -maxhits 1 -uc ${DATASET}_unique.uc -mothur_shared_out ${DATASET}_uniqueOtuTab.shared -otutabout ${DATASET}_uniqueOtuTab.txt

cp ${DATASET}_uniqueOtuTab.* ${RESULTSDIR}
echo -e `grep -v -c "*" ${DATASET}_unique.uc`" merged unique read pairs map to clusters" >> ${RUNDIR}${DATASET}.log
echo -e "-- \n" >> ${RUNDIR}${DATASET}.log
# echo -e "-- \n" >> ${RUNDIR}${DATASET}.log

# Classify sequences 
echo -e "Classify sequences:" >> ${RUNDIR}${DATASET}.log
echo -e "------------------" >> ${RUNDIR}${DATASET}.log
if [[ $CLASSIFY == "Y" ]]
then
 echo -e "Classifying OTUs:" >> ${RUNDIR}${DATASET}.log
 cp ${RESULTSDIR}${DATASET}_unique.fa ${DATASET}_unique.fa
 ${MOTHUR}mothur "#classify.seqs(fasta=${DATASET}_unique.fa,template=silva.nr_v119.align, taxonomy=silva.nr_v119.tax, processors=10)" > tmp.log
 grep -m 1 'mothur >' tmp.log >> ${RUNDIR}${DATASET}.log 
 grep -m 1 -A1000 'Output File Names:' tmp.log | head -n -3  >> ${RUNDIR}${DATASET}.log   
 echo -e "-- \n" >> ${RUNDIR}${DATASET}.log
 cp ${DATASET}_unique.nr_v119.wang.taxonomy ${RESULTSDIR}${DATASET}_unique_silva.nrv119.taxonomy
fi
cp 

# Return results 
cp ${RUNDIR}${DATASET}.log ${RESULTSDIR}
tar -zcvf ${RUNDIR}${DATASET}.OTU.tar.gz ${RESULTSDIR}/*
#rm *.logfile silva.nr* tmp.log
cd ..
#tar -zcvf$ ${RUNDIR}${DATASET}.OTU.tar.gz ${WORKDIR}/*

# Cleanup
rm -r ${WORKDIR}
rm -r ${RESULTSDIR}
