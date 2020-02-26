#!/bin/bash

AUTHOR="imandric1"



################################################################
##########          The main template script          ##########
################################################################

toolName="fiona"
toolPath="/u/home/n/ngcrawfo/project-zarlab/igor/VIRAL_CORRECTION/tools/fiona-0.2.8-Linux-x86_64/bin/fiona"

# STEPS OF THE SCRIPT
# 1) prepare input if necessary
# 2) run the tool
# 3) transform output if necessary
# 4) compress output


# THE COMMAND LINE INTERFACE OF THE WRAPPER SCRIPT
# $tool $input1 $input2 $outdir $kmers $others
# |      mandatory part       | | extra part |
# <---------------------------> <------------>




if [ $# -lt 4 ]
then
echo "********************************************************************"
echo "Script was written for project : Best practices for conducting benchmarking in the most comprehensive and reproducible way"
echo "This script was written by Igor Mandric"
echo "********************************************************************"
echo ""
echo "1 <input1> - _1.fastq"
echo "2 <input2> - _2.fastq"
echo "3 <outdir> - dir to save the output"
echo "4 <glen>   - genome length"
echo "--------------------------------------"
exit 1
fi



# mandatory part
input1=$1
input2=$2
outdir=$3

# extra part (tool specific)
glen=$4



# STEP 0 - create output directory if it does not exist

mkdir -p $outdir
pwd=$PWD
cd $outdir
outdir=$PWD
cd $pwd
logfile=$outdir/report_$(basename ${input1%.*})_${toolName}_${kmer}.log

# -----------------------------------------------------

echo "START" >> $logfile

# STEP 1 - prepare input if necessary (ATTENTION: TOOL SPECIFIC PART!)

cat $input1 $input2 > $outdir/merged_input_file.fastq

# -----------------------------------





# STEP 2 - run the tool (ATTENTION: TOOL SPECIFIC PART!)

now="$(date)"
printf "%s --- RUNNING %s\n" "$now" $toolName >> $logfile

# run the command
res1=$(date +%s.%N)
$toolPath -g $glen $outdir/merged_input_file.fastq $outdir/one_output_file.fasta >> $logfile 2>&1
res2=$(date +%s.%N)
dt=$(echo "$res2 - $res1" | bc)
dd=$(echo "$dt/86400" | bc)
dt2=$(echo "$dt-86400*$dd" | bc)
dh=$(echo "$dt2/3600" | bc)
dt3=$(echo "$dt2-3600*$dh" | bc)
dm=$(echo "$dt3/60" | bc)
ds=$(echo "$dt3-60*$dm" | bc)
now="$(date)"
printf "%s --- TOTAL RUNTIME: %d:%02d:%02d:%02.4f\n" "$now" $dd $dh $dm $ds >> $logfile

now="$(date)"
printf "%s --- FINISHED RUNNING %s %s\n" "$now" $toolName >> $logfile

# ---------------------




# STEP 3 - transform output if necessary (ATTENTION: TOOL SPECIFIC PART!)

# if you need to transform fasta to fastq - here is the command:
#     awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' file.fastq > file.fasta

now="$(date)"
printf "%s --- TRANSFORMING OUTPUT\n" "$now" >> $logfile

awk 'NR%2{printf "%s ",$0;next;}1' $outdir/one_output_file.fasta | awk '{lastfield=$NF;$NF="";print "@"substr($0,2,length($0)-1)"\n"lastfield"\n+\n"lastfield}' | gzip > $outdir/${toolName}_$(basename ${input1%.*})_${glen}.corrected.fastq.gz

now="$(date)"
printf "%s --- TRANSFORMING OUTPUT DONE\n" "$now" >> $logfile


# remove intermediate files
rm $outdir/merged_input_file.fastq
rm $outdir/one_output_file.*

# --------------------------------------



printf "DONE" >> $logfile



