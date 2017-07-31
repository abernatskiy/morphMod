#!/bin/bash

# A scripts to generate the final scripts for the gigantic convergence experiment

generateOneScript() {
	local base_script="${1}"
	local num_segments="${2}"
	local stop_after_generations="${3}"
	local pareto_logging_period="${4}"
	local queue="${5}"
	local requested_time="${6}"
	local max_jobs="${7}"
	local simulation_time="${8}"
	local simulation_time_step="${9}"
	local stages_multiplier="${10}"
	local indiv_classes="${11}"
	cat "${base_script}.py" | sed \
		-e "s/segments = [0-9]*/segments = ${num_segments}/" \
		-e "s/'genStopAfter': [0-9]*/'genStopAfter': ${stop_after_generations}/" \
		-e "s/'logParetoFrontPeriod': [0-9]*/'logParetoFrontPeriod': ${pareto_logging_period}/" \
		-e "s/queue = '[a-z]*'/queue = '${queue}'/" \
		-e "s/maxJobs = [0-9]*/maxJobs = ${max_jobs}/" \
		-e "s/expectedWallClockTime = '.*'/expectedWallClockTime = '${requested_time}'/" \
		-e "s/'simulationTime': [\.0-9]*/'simulationTime': ${simulation_time}/" \
		-e "s/'timeStep': [\.0-9]*/'timeStep': ${simulation_time_step}/" \
		-e "s/mult = [0-9]*/mult = ${stages_multiplier}/" \
		-e "s/'compositeClass0', \[.*\]/'compositeClass0', \[${indiv_classes}\]/" > "autoscripts/${base_script}_N${num_segments}.py"
}

generateFromOneBase() {
	local base_script="${1}"
	# Regarding convergence times: about 100 gens for 3 segments at probability of morphological mutation of 0.2, 400 for 5, 10 gets stuck indefinitely.
	# The times are doubled where known and for 10 segments a placeholder 2000 is used for now
	generateOneScript "${base_script}" 3 125 1 'shortq' '03:00:00' 2 '3.0' '0.05' 5 "'integerVectorSymmetricRangeMutations', 'integerVectorRandomJumps'"
	generateOneScript "${base_script}" 5 625 5 'shortq' '03:00:00' 8 '3.0' '0.05' 5 "'integerVectorSymmetricRangeMutations', 'integerVectorRandomJumps'"
	generateOneScript "${base_script}" 10 1000 10 'workq' '30:00:00' 100 '6.0' '0.1' 10 "'integerVectorSymmetricRangeMutations'"
}

generateFromOneBase simpleTimeSeries
generateFromOneBase rateSwipe
generateFromOneBase rateSizeSwipe
generateFromOneBase ccSwipe
