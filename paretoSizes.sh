#!/bin/bash

for dir in */; do
	for file in ./${dir}/*/paretoFront_gen*.log ; do
		wc -l $file
	done | cut -d ' ' -f1 | sort -n | uniq -c > ${dir}/paretoHistogram
done
