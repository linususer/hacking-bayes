SHELL := /bin/bash
prepareSimulations:
	mkdir -p data
	R -f createDB.R

makeAllSimulations:
	@echo "This will take a long time. Consider simulating with less parameters or only specific simulations."
	@echo "If data/results.db exists, then all data will be overwritten."
	@read -p "Do you want to continue? [Y/N]: " confirm; \
	if [ "$$confirm" = "n" ] || [ "$$confirm" = "N" ]; then exit 1; \
	elif [ "$$confirm" = "y" ] || [ "$$confirm" = "Y" ]; then exit 0; \
	fi
	R -f calculations/application-sim.R
	@echo "Done with application example simulations.\n"
	@echo "Start cauchy prior simulation.\n"
	R -f calculations/cauchy-functions.R
	@echo "Done with cauchy prior simulation.\n"
	@echo "Start n_start parameter simulation with normal prior.\n"
	R -f calculations/rouder-simulations.R
	@echo "Done with n_start parameter simulation with normal prior.\n"
	@echo "Start Realistic simulation with cauchy prior\n"
	R -f calculations/realistic-sim-par.R
	@echo "Done with Realistic simulation with cauchy prior\n"

exportSimulations:
	@echo "Change name of database\n"
	mv data/results.duckdb data/hacking-bayes-factors.duckdb
	@echo "Create an aggregated summary"