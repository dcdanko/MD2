from unittest import TestCase

from microbe_directory import annotate_taxa


class TestTaxaAnnotate(TestCase):

    def test_annotate(self):
        taxa = [
            'Cutibacterium acnes', 'Cutibacterium granulosum',
            'Bacteria', 'Escherichia', 'Escherichia coli',
            'Cutibacterium'
        ]
        tbl = annotate_taxa(taxa)
        self.assertEqual(len(taxa), tbl.shape[0])
        self.assertGreater(tbl.shape[1], 0)
