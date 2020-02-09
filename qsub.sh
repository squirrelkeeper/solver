#!/bin/bash 

timestamp="q"`date +%y-%m-%d`
sub_timestamp=`date +%H-%M-%S`
folder_path="data/"$timestamp
sub_folder_path=$folder_path"/"$sub_timestamp

argfile="$sub_folder_path/argfile_$timestamp-$sub_timestamp.arg"
jobfile="jobfile_$timestamp-$sub_timestamp.out"

prog="mll_oefb_v3"
inc="inc"

par_str="-option1 -par1 500 -par2"

sim_start=0
sim_end=8
sim_num=1

mem=1
prio=1

echo timestamp: $timestamp

mkdir -p $sub_folder_path

for (( i=0; i<$sim_num; i++ ))
do
	current_par=$(echo "scale=8;$i*($sim_end-$sim_start)/$sim_num + sim_start" | bc -l | awk '{printf "%f", $0}')
	echo "$par_str $current_par" >> $argfile
done

qsub -mem $mem -p $prio -m n -w $sub_folder_path -o $jobfile -argfile $argfile $prog

zip $sub_folder_path"/"backup.zip -r $prog $prog.cpp qsub.sh $inc makefile $argfile
