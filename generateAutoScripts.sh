#!/bin/bash

# A scripts to generate the final scripts for the gigantic convergence experiment

generateOneScript() {
	local base_script="${1}"
	local num_segments="${2}"
	local stop_after_generations="${3}"
	local pareto_logging_period="${4}"
	local queue="${5}"
	local requested_time="${6}"
	cat "${base_script}.py" | sed \
		-e "s/segments = [0-9]*/segments = ${num_segments}/" \
		-e "s/'genStopAfter': [0-9]*/'genStopAfter': ${stop_after_generations}/" \
		-e "s/'logParetoFrontPeriod': [0-9]*/'logParetoFrontPeriod': ${pareto_logging_period}/" \
		-e "s/queue = '[a-z]*'/queue = '${queue}'/" \
		-e "s/expectedWallClockTime = '.*'/expectedWallClockTime = '${requested_time}'/" > "autoscripts/${base_script}_N${num_segments}.py"
}

generateFromOneBase() {
	local base_script="${1}"
	# Regarding convergence times: about 100 gens for 3 segments at probability of morphological mutation of 0.2, 400 for 5, 10 gets stuck indefinitely.
	# The times are doubled where known and for 10 segments a placeholder 2000 is used for now
	generateOneScript "${base_script}" 3 125 1 'shortq' '03:00:00'
	generateOneScript "${base_script}" 5 625 5 'shortq' '03:00:00'
	generateOneScript "${base_script}" 10 2000 10 'workq' '30:00:00'
}

generateFromOneBase simpleTimeSeries
generateFromOneBase rateSwipe
generateFromOneBase rateSizeSwipe
