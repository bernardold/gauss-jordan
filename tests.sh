#!/bin/bash
log_all="all.csv"
echo "proc_num, dimension, groups, time" >> ${log_all}
for p in 2 4 8 16
do
	for k in 4 16 32 64
	do
		let dimension=p*k
		echo "Start test - proc_num[" ${p} "] and dimension [" ${dimension} "]"
		for groups in 0 1
		do
			echo -n "${p}, ${dimension}, ${groups}, ">> ${log_all}
			echo `mpiexec -n ${p} ./gj ${dimension} ${groups}` >> ${log_all}
		done
	done
done
echo "Done!"
