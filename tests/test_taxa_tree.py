from unittest import TestCase

from microbe_directory.taxa_tree import NCBITaxaTree


class TestTaxaTree(TestCase):

    def setUp(self):
        """Test that we can make a taxa tree."""
        self.tree = NCBITaxaTree.parse_files()

    def test_tree_sort(self):
        taxa = [
            'Cutibacterium acnes', 'Cutibacterium granulosum',
            'Bacteria', 'Escherichia', 'Escherichia coli',
            'Cutibacterium'
        ]
        sort = self.tree.taxa_sort(taxa)
        self.assertEqual(len(taxa), len(sort))

    def test_get_phylum(self):
        phylum = self.tree.phylum('Escherichia coli')
        self.assertEqual(phylum, 'Proteobacteria')

    def test_get_rank(self):
        rank = self.tree.rank('Escherichia coli')
        self.assertEqual(rank, 'species')

    def test_get_parent(self):
        parent = self.tree.parent('Cutibacterium acnes')
        self.assertEqual(parent, 'Cutibacterium')

    def test_get_ancestors(self):
        ancestors = self.tree.ancestors('Cutibacterium acnes')
        true = [
            'Cutibacterium acnes',
            'Cutibacterium',
            'Propionibacteriaceae',
            'Propionibacteriales',
            'Actinobacteria',
            'Actinobacteria',
            'Terrabacteria group',
            'Bacteria',
            'cellular organisms',
            'root'
        ]
        for i, anc in enumerate(ancestors):
            self.assertEqual(anc, true[i])

    def test_place_virus(self):
        taxon_id, rank, taxon_superkingdom = self.tree.place_microbe('Flaviviridae')
        self.assertEqual(taxon_id, '11050')
        self.assertEqual(rank, 'family')
        self.assertEqual(taxon_superkingdom, 'Viruses')

    def test_place_bacteria(self):
        taxon_id, rank, taxon_superkingdom = self.tree.place_microbe('Escherichia coli')
        self.assertEqual(taxon_id, '562')
        self.assertEqual(rank, 'species')
        self.assertEqual(taxon_superkingdom, 'Bacteria')

    def test_place_euk(self):
        taxon_id, rank, taxon_superkingdom = self.tree.place_microbe('Saccharomyces cerevisiae')
        self.assertEqual(taxon_id, '4932')
        self.assertEqual(rank, 'species')
        self.assertEqual(taxon_superkingdom, 'Eukaryota')
