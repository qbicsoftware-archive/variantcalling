#!/bin/sh
#PBS -q cfc
#PBS -A qbic
#PBS -l nodes=1:ppn=10:cfc
#PBS -l walltime=200:00:00
#PBS -l mem=50g
# properties = {properties}

set -e

echo Running on machine $(hostname)

module load qbic/anaconda
module load devel/java_jdk/1.7.0u45
module load qbic/samtools/1.3
module load qbic/ngs-bits
module load qbic/bwa
module load qbic/stampy
module load qbic/picard/git
module load bio/gatk/3.3
module load qbic/freebayes/0.9
module load qbic/vcflib
module load qbic/bcftools
module load qbic/vcftools

{exec_job}
exit 0
