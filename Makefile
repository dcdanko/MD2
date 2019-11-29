.DEFAULT_GOAL: help

help:
	@echo "clean - remove all build, test, coverage and Python artifacts"
	@echo "virus - build virus table"
	@echo "bact - build bacteria table"
	@echo "fungi - build fungi table"
	@echo "all - build bacteria, virus, and fungi tables"
	@echo "test - run unit tests"

clean: 
	rm -r results | exit 0

virus:
	mkdir -p results/virus
	microbe_directory taxa-table -s viruses > results/virus/viruses_taxa.csv
	microbe_directory merge-csvs results/virus/viruses_taxa.csv datasets/virus/*.csv > results/virus/viruses_merged.csv
	microbe_directory clean-file results/virus/viruses_merged.csv > results/virus/viruses_cleaned.csv
	cp results/virus/viruses_cleaned.csv results/viruses.csv

bact:
	mkdir -p results/bacts
	microbe_directory taxa-table -s bacteria > results/bacts/bacteria_taxa.csv
	microbe_directory merge-csvs results/bacts/bacteria_taxa.csv datasets/bact/*.csv > results/bacts/bacteria_merged.csv
	microbe_directory clean-file results/bacts/bacteria_merged.csv > results/bacts/bacteria_cleaned.csv
	cat results/bacts/bacteria_cleaned.csv | grep -Fv 'sp.' | grep -iFv 'culture' | grep -iFv 'candidate' | grep -iFv 'candidatus' > results/bacts/bacteria_filtered.csv
	microbe_directory infill-bacteria results/bacts/bacteria_filtered.csv > results/bacts/bacteria_filled.csv
	cp results/bacts/bacteria_filled.csv results/bacteria.csv

euks:
	mkdir -p results/euks
	microbe_directory taxa-table -s eukaryota > results/euks/eukaryota_taxa.csv
	microbe_directory merge-csvs results/euks/eukaryota_taxa.csv datasets/euk/*.csv > results/euks/eukaryota_merged.csv
	microbe_directory clean-file results/euks/eukaryota_merged.csv > results/euks/eukaryota_cleaned.csv
	cp results/euks/eukaryota_cleaned.csv results/eukaryota.csv

all: virus euks bact

test:
	pytest tests
