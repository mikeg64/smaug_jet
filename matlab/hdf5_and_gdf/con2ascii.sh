#!/bin/bash

#$ -l mem=16G
#$ -l rmem=16G
#$ -N con2h5
#$ -j y
##$ -t 5-8:1
#$ -V 

#$ -l h_rt=8:00:00


##module add libs/cuda/4.0.17
#module add libs/cuda/6.5.14
module load apps/matlab/2018b


#export TIMECOUNTER=0
#source timeused
   matlab -nosplash -r "sac3Dtoascii"
#source timeused

#cd /fastdata/cs1mkg/smaug_wash/washmc_2p5_2p5_12p5_mach180_uni6

#tar -zcvf h5.tgz h5
