.DEFAULT_GOAL: help

help:
	@echo "clean - remove all build, test, coverage and Python artifacts"
	@echo "virus - build virus table"
	@echo "bact - build bacteria table"
	@echo "fungi - build fungi table"
	@echo "all - build bacteria, virus, and fungi tables"
	@echo "test - run unit tests"

clean: 
	rm -r build

virus:
	mkdir -p build/virus
	microbe_directory taxa-table -s viruses > build/virus/viruses_taxa.csv
	microbe_directory merge-csvs build/virus/viruses_taxa.csv datasets/*.csv > build/virus/viruses_merged.csv
	microbe_directory clean-file build/virus/viruses_merged.csv > build/virus/viruses_cleaned.csv

bact:
	mkdir -p build/bacts
	microbe_directory taxa-table -s bacteria > build/bacts/bacteria_taxa.csv
	microbe_directory merge-csvs build/bacts/bacteria_taxa.csv datasets/*.csv > build/bacts/bacteria_merged.csv
	microbe_directory infill-bacteria build/bacts/bacteria_merged.csv > build/bacts/bacteria_filled.csv
	microbe_directory clean-file build/bacts/bacteria_filled.csv > build/bacts/bacteria_cleaned.csv

euks:
	mkdir -p build/euks
	microbe_directory taxa-table -s eukaryota > build/euks/eukaryota_taxa.csv
	microbe_directory merge-csvs build/euks/eukaryota_taxa.csv datasets/*.csv > build/euks/eukaryota_merged.csv
	microbe_directory clean-file build/euks/eukaryota_merged.csv > build/euks/eukaryota_cleaned.csv

all: virus fungi bact

test:
	pytest tests
