
#PBS -o results/output.txt
#PBS -e results/errors.txt 
#PBS -l walltime=8:0:0
#PBS -l nodes=1:ppn=2
#PBS -q intel
#PBS -l mem=4gb
. 
mpirun -np 2 python results/microcircuit.py
