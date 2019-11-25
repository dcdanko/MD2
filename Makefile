.DEFAULT_GOAL: help

help:
	@echo "clean - remove all build, test, coverage and Python artifacts"
	@echo "virus - build virus table"
	@echo "bact - build bacteria table"
	@echo "fungi - build fungi table"
	@echo "all - build bacteria, virus, and fungi tables"
	@echo "test - run unit tests"

clean: 
	rm -r results

virus:
	mkdir -p results/virus
	microbe_directory taxa-table -s viruses > results/virus/viruses_taxa.csv
	microbe_directory merge-csvs results/virus/viruses_taxa.csv datasets/*.csv > results/virus/viruses_merged.csv
	microbe_directory clean-file results/virus/viruses_merged.csv > results/virus/viruses_cleaned.csv

bact:
	mkdir -p results/bacts
	microbe_directory taxa-table -s bacteria > results/bacts/bacteria_taxa.csv
	microbe_directory merge-csvs results/bacts/bacteria_taxa.csv datasets/*.csv > results/bacts/bacteria_merged.csv
	microbe_directory infill-bacteria results/bacts/bacteria_merged.csv > results/bacts/bacteria_filled.csv
	microbe_directory clean-file results/bacts/bacteria_filled.csv > results/bacts/bacteria_cleaned.csv

euks:
	mkdir -p results/euks
	microbe_directory taxa-table -s eukaryota > results/euks/eukaryota_taxa.csv
	microbe_directory merge-csvs results/euks/eukaryota_taxa.csv datasets/*.csv > results/euks/eukaryota_merged.csv
	microbe_directory clean-file results/euks/eukaryota_merged.csv > results/euks/eukaryota_cleaned.csv

all: virus fungi bact

test:
	pytest tests
