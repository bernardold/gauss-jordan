#!/bin/bash
log_all="all.csv"
echo "proc_num, dimension, groups, time" >> ${log_all}
for groups in 1 0
do
	for k in 64 128 256 512
	do
		for p in 2 4 8 16 32
		do
			let dimension=p*k
			echo "Start test - proc_num[ ${p} ], dimension [ ${dimension} ], use k-groups[ ${groups} ] at `date +%d/%m-%R`"
			echo -n "${p}, ${dimension}, ${groups}, ">> ${log_all}
			echo `mpiexec -n ${p} ./gj ${dimension} ${groups}` >> ${log_all}
		done
	done
done
echo "Done!"
