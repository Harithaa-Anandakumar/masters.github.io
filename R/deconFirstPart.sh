#!/bin/bash
#SBATCH --job-name="decon7"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mail-user=harithaa.anand.125@gmail.com
#SBATCH --time=0-48:00
#SBATCH -o decon7.o%j
#SBATCH -e decon7.e%j


/usr/users/hananda/workspace/thesis/decontaMiner_1.4/shell_scripts/decontaMiner.sh \
 -i /usr/users/hananda/rawdata/RnaSeq/thesisBam/unmapped/seven \
 -o /usr/users/hananda/rawdata/RnaSeq/thesisBam/decon7 \
 -c /usr/users/hananda/workspace/copied_decon/decontaMiner_1.4/config_files/rna/configure.txt \
 -s S -bv -Q n