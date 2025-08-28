#!/bin/bash
#SBATCH --job-name=BL3_Simulation    ### Job Name
#SBATCH --output=main.out       ### File in which to store job output
#SBATCH --error=main.err        ### File in which to store job error messages
#SBATCH --partition=centos7      ### Partition (default is 'defq')
#SBATCH --qos=normal          ### Quality of Service (like a queue in PBS)
#SBATCH --time=0-00:010:00     ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1             ### Node count required for the job
#SBATCH --ntasks-per-node=1   ### Number of tasks to be launched per Node

module load gcc/9.5.0
g++ -O2 -std=c++17 main.c -o main
./main