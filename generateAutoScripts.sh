#!/bin/bash

# A scripts to generate the final scripts for the gigantic convergence experiment

generateOneScript() {
	local base_script="${1}"
	local num_segments="${2}"
	local stop_after_generations="${3}"
	cat "${base_script}.py" | sed -e "s/segments = [0-9]*/segments = ${num_segments}/" | sed -e "s/'genStopAfter': [0-9]*/'genStopAfter': ${stop_after_generations}/" > "autoscripts/${base_script}_${num_segments}segments.py"
}

generateFromOneBase() {
	local base_script="${1}"
	generateOneScript "${base_script}" 3 200
	generateOneScript "${base_script}" 5 800
	generateOneScript "${base_script}" 10 1000
}

generateFromOneBase simpleTimeSeries
