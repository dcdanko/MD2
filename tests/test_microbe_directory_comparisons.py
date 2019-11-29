"""Test suite for experimental functions."""

import pandas as pd

from unittest import TestCase
from os.path import join, dirname

from microbe_directory.comparisons import (
    compare_microbe_directory_dataframes,
    compare_taxa_lists,
    compare_taxa_lists_abundances
)
from microbe_directory.comparisons.statistics import (
    compare_categorical,
    compare_numeric,
    compare_numeric_abundances,
    compare_categorical_abundances,
)
from microbe_directory.comparisons.constants import (
    MICROBE_DIRECTORY,
)

PACKET_DIR = join(dirname(__file__), 'built_packet')


class TestMicrobeDirectoryComparisons(TestCase):
    """Test suite for comparing taxa lists using the microbe directory."""

    def test_compare_categorical_binary(self):
        """Test that we can run UMAP."""
        categorical_test = compare_categorical(
            'yes',
            pd.Series(['yes', 'yes', 'no', 'yes', 'yes']),
            pd.Series(['no', 'no', 'no', 'yes', 'yes', 'no', 'no']),
        )
        self.assertTrue(0 <= categorical_test['p-value'].all() <= 1)
        self.assertTrue(0 < categorical_test['abundance_in'].all())
        self.assertTrue(0 < categorical_test['abundance_out'].all())

    def test_compare_categorical_multi(self):
        """Test that we can run UMAP."""
        categorical_test = compare_categorical(
            'A',
            pd.Series(['A', 'B', 'D']),
            pd.Series(['B', 'C', 'D', 'E'])
        )
        self.assertTrue(0 <= categorical_test['p-value'].all() <= 1)
        self.assertTrue(0 < categorical_test['abundance_in'].all())
        self.assertTrue(0 < categorical_test['abundance_out'].all())

    def test_compare_numeric(self):
        """Test that we can run fractal."""
        numeric_test = compare_numeric(
            pd.Series([0, 1, 3, 0, 1, 1, 2, 2]),
            pd.Series([2, 2, 1, 3, 1, 3, 4]),
        )
        self.assertTrue(0 <= numeric_test['p-value'] <= 1)
        self.assertTrue(0 < numeric_test['abundance_in'])
        self.assertTrue(0 < numeric_test['abundance_out'])

    def test_compare_dataframes(self):
        dataframe_test = compare_microbe_directory_dataframes(
            pd.DataFrame(MICROBE_DIRECTORY.iloc[0:5, 7:30]),
            pd.DataFrame(MICROBE_DIRECTORY.iloc[9:14, 7:30]),
        )
        self.assertTrue(len(dataframe_test.columns) == 7)

    def test_compare_taxa_lists(self):
        taxa_list_test = compare_taxa_lists(
            MICROBE_DIRECTORY.iloc[0:5].index.tolist(),
            MICROBE_DIRECTORY.iloc[9:14].index.tolist(),
        )
        self.assertTrue(len(taxa_list_test.columns) == 7)

    def test_compare_categorical_abundances(self):
        cat_abundances_test = compare_categorical_abundances(
            'A',
            {'A': 0.2, 'B': 0.3, 'C': 0.5},
            {'B': 0.25, 'C': 0.4, 'D': 0.25, 'E': 0.1},
        )
        self.assertTrue(0 <= cat_abundances_test['p-value'] <= 1)
        self.assertTrue(0 < cat_abundances_test['abundance_in'].all())
        self.assertTrue(0 < cat_abundances_test['abundance_out'].all())

    def test_compare_numeric_abundances(self):
        numeric_abundances_test = compare_numeric_abundances(
            {5: 0.2, 6: 0.25, 7: 0.25, 8: 0.1, 9: 0.2},
            {4: 0.2, 6: 0.125, 7: 0.3, 8: 0.375},
        )
        self.assertTrue(0 <= numeric_abundances_test['p-value'] <= 1)
        self.assertTrue(0 < numeric_abundances_test['abundance_in'])
        self.assertTrue(0 < numeric_abundances_test['abundance_out'])

    def test_compare_taxa_lists_abundances(self):
        taxa_list_abundances_test = compare_taxa_lists_abundances(
            pd.Series({
                'Bacteriovorax marinus': 0.2, 'Staphylococcus phage 44AHJD': 0.1,
                'Junonia coenia densovirus': 0.3, 'Listeria phage B025': 0.25,
                'Prevotella copri': 0.15}
            ),
            pd.Series({
                'Prevotella nigrescens': 0.25, 'Mycobacterium phage PattyP': 0.3,
                'Clostridium asparagiforme': 0.15, 'Nocardiopsis halophila': 0.2,
                'Cucumber Bulgarian virus': 0.1
            })
        )
        self.assertTrue(len(taxa_list_abundances_test.columns) == 7)
