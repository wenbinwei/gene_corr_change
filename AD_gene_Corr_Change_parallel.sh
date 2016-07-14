#!/bin/bash
# To run this script on Iceberg using the following 
# cd /shared/hidelab2/user/md1wwxx/ADCorrChange/src/
# qsub AD_gene_Corr_Change_parallel.sh
#$ -pe openmp 3
#memory requests are per-core
#$ -l rmem=14G -l mem=14G
#Prefer the hidelab queue but spill over to over queues if it is full
#$ -P hidelab
#$ -m bea # send mails at beginning, end and if aborted unexpectedly
#$ -M W.Wei@Sheffield.ac.uk # mail sent to this address
#$ -o AD_gene_Corr_Change_parallel.log # output file 
#$ -j y # send output and error to same file

cd /shared/hidelab2/user/md1wwxx/ADCorrChange/src/
# load R
module load apps/R
R CMD BATCH AD_gene_Corr_Change_parallel.R AD_gene_Corr_Change_parallel.Rout