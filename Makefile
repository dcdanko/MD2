.DEFAULT_GOAL: help

help:
	@echo "clean - remove all build, test, coverage and Python artifacts"
	@echo "virus - build virus table"
	@echo "bact - build bacteria table"
	@echo "euks - build fungi table"
	@echo "all - build bacteria, virus, and fungi tables"
	@echo "metasub - preprocess metasub tables"
	@echo "store - copy tables in results to source dir for distribution"
	@echo "test - run unit tests"

clean: 
	rm -r results | exit 0

virus:
	mkdir -p results/virus
	microbe_directory taxa-table -s viruses > results/virus/viruses_taxa.csv
	microbe_directory merge-csvs results/virus/viruses_taxa.csv datasets/virus/*.csv datasets/shared/*.csv \
		> results/virus/viruses_merged.csv
	microbe_directory clean-file results/virus/viruses_merged.csv > results/virus/viruses_cleaned.csv
	cp results/virus/viruses_cleaned.csv results/viruses.csv

bact:
	mkdir -p results/bacts
	microbe_directory taxa-table -s bacteria > results/bacts/bacteria_taxa.csv
	microbe_directory merge-csvs results/bacts/bacteria_taxa.csv datasets/bact/*.csv datasets/shared/*.csv \
		> results/bacts/bacteria_merged.csv
	microbe_directory clean-file results/bacts/bacteria_merged.csv > results/bacts/bacteria_cleaned.csv
	cat results/bacts/bacteria_cleaned.csv \
		| grep -Fv 'sp.' \
		| grep -iFv 'culture' \
		| grep -iFv 'candidate' \
		| grep -iFv 'candidatus' \
		> results/bacts/bacteria_filtered.csv
	microbe_directory infill-bacteria results/bacts/bacteria_filtered.csv > results/bacts/bacteria_filled.csv
	microbe_directory composite-bacteria results/bacts/bacteria_filled.csv > results/bacts/bacteria_composite.csv
	cp results/bacts/bacteria_composite.csv results/bacteria.csv

euks:
	mkdir -p results/euks
	microbe_directory taxa-table -s eukaryota > results/euks/eukaryota_taxa.csv
	microbe_directory merge-csvs results/euks/eukaryota_taxa.csv datasets/euk/*.csv datasets/shared/*.csv \
		> results/euks/eukaryota_merged.csv
	microbe_directory clean-file results/euks/eukaryota_merged.csv > results/euks/eukaryota_cleaned.csv
	cp results/euks/eukaryota_cleaned.csv results/eukaryota.csv

metasub:
	microbe_directory preprocess metasub -f city \
		datasets/needs_preprocessing/metasub/refseq.krakenhll_species.csv \
		datasets/needs_preprocessing/metasub/complete_metadata.csv \
		> datasets/shared/metasub_city.csv
	microbe_directory preprocess metasub -f continent \
		datasets/needs_preprocessing/metasub/refseq.krakenhll_species.csv \
		datasets/needs_preprocessing/metasub/complete_metadata.csv \
		> datasets/shared/metasub_continent.csv
	microbe_directory preprocess metasub -f surface_ontology_fine \
		datasets/needs_preprocessing/metasub/refseq.krakenhll_species.csv \
		datasets/needs_preprocessing/metasub/complete_metadata.csv \
		> datasets/shared/metasub_surface.csv

store:
	gzip -c results/eukaryota.csv > microbe_directory/stored_final_tables/eukaryota.csv.gz
	gzip -c results/bacteria.csv > microbe_directory/stored_final_tables/bacteria.csv.gz
	gzip -c results/viruses.csv > microbe_directory/stored_final_tables/viruses.csv.gz


all: virus euks bact

test:
	pytest tests
