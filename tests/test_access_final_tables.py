from unittest import TestCase

from microbe_directory import (
    bacteria,
    eukaryote,
    virus,
    md1,
)


class TestAccessFinalTables(TestCase):

    def test_get_virus(self):
        tbl = virus()
        self.assertGreater(tbl.shape[0], 0)
        self.assertGreater(tbl.shape[1], 0)

    def test_get_eukaryote(self):
        tbl = eukaryote()
        self.assertGreater(tbl.shape[0], 0)
        self.assertGreater(tbl.shape[1], 0)

    def test_get_bacteria(self):
        tbl = bacteria()
        self.assertGreater(tbl.shape[0], 0)
        self.assertGreater(tbl.shape[1], 0)

    def test_get_md1(self):
        tbl = md1()
        self.assertGreater(tbl.shape[0], 0)
        self.assertGreater(tbl.shape[1], 0)
