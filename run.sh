#!/bin/bash

#Author: Siqi Wang
#File: run.sh
#Create Date: 2017-11-23 11:00

snakemake -p --rerun-incomplete -k -j 128 --cluster "bsub -n {cluster.ppn} -J {params.jobname} -q {cluster.queue} -o log/{params.jobname}.pbs.out -e log/{params.jobname}.pbs.err -R \"span[hosts=1]\"" --jobscript jobscript.pbs --jobname "{rulename}.{jobid}.pbs" --cluster-config cluster.json 2>run.log
