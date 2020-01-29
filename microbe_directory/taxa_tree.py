from os import environ
from os.path import join, dirname
import gzip
import glob
import pandas as pd
from .dataset_modification import (
    taxa_to_organism,
)
from .constants import (
    RANK_LIST,
    ALLOWED_SUPERKINGDOMS,
    CANONICAL_RANKS,
    ROOT_RANK,
)


NCBI_DELIM = '\t|'  # really...
NAMES_ENV_VAR = 'MD2_NCBI_NAMES'
NODES_ENV_VAR = 'MD2_NCBI_NODES'
RANKEDLINEAGE_ENV_VAR = 'MD2_NCBI_RANKEDLINEAGE'
NAMES_DEF = join(dirname(__file__), 'ncbi_tree/names.dmp.gz')
NODES_DEF = join(dirname(__file__), 'ncbi_tree/nodes.dmp.gz')
RANKEDLINEAGE_DEF = join(dirname(__file__), 'ncbi_tree/rankedlineage.dmp.gz')


class TaxonomicRankError(Exception):
    pass


class NCBITaxaTree:

    def __init__(self, parent_map, names_to_nodes, nodes_to_name):
        self.parent_map = parent_map
        self.names_to_nodes = names_to_nodes
        self.nodes_to_name = nodes_to_name

    def _node(self, taxon_name):
        return self.names_to_nodes[taxon_name]

    def _name(self, node_num):
        return self.nodes_to_name[node_num]['name']

    def rank(self, taxon_name):
        return self.nodes_to_name[self._node(taxon_name)]['rank']

    def taxon_id(self, taxon_name):
        return self._node(taxon_name)

    def parent(self, taxon):
        """Return the name of the parent taxon."""
        return self._name(self.parent_map[self._node(taxon)])

    def ancestor_rank(self, rank, taxon, default=None):
        """Return the ancestor of taxon at the given rank."""
        try:
            parent_num = self.parent_map[self._node(taxon)]
            while int(parent_num) > 1:
                if rank == self.nodes_to_name[parent_num]['rank']:
                    return self.nodes_to_name[parent_num]['name']
                parent_num = self.parent_map[parent_num]
        except KeyError:
            if default is None:
                raise
        return default

    def ancestors(self, taxon, max_rank=ROOT_RANK):
        return self.ancestors_list(taxon, max_rank=max_rank)

    def ancestors_list(self, taxon, max_rank=ROOT_RANK):
        """Return a phylogenetically sorted list of ancestors of taxon including taxon."""
        max_rank_index = RANK_LIST.index(max_rank)
        try:
            rank = self.rank(taxon)
            taxon_rank_index = RANK_LIST.index(rank)
        except ValueError:
            raise TaxonomicRankError(f'Requested rank {rank} is not in rank list.')
        if taxon_rank_index > max_rank_index:
            raise TaxonomicRankError(f'Requested rank {rank} is above {taxon}.')
        parent_num = self.parent_map[self._node(taxon)]
        parent_rank = self.nodes_to_name[parent_num]['rank']
        try:
            rank_index = RANK_LIST.index(parent_rank)
        except ValueError:
            rank_index = -1
        ancestor_name_list = [taxon]
        while max_rank_index > rank_index:
            ancestor_name_list.append(self.nodes_to_name[parent_num]['name'])
            if int(parent_num) == 1:  # root
                break
            parent_num = self.parent_map[parent_num]
            parent_rank = self.nodes_to_name[parent_num]['rank']
            try:
                rank_index = RANK_LIST.index(parent_rank)
            except ValueError:
                rank_index = -1
        return ancestor_name_list

    def canonical_taxonomy(self, taxon):
        """Return a dict with the canonical (KPCOFGS) taxonomy for a taxon and the taxon id."""
        parent_num = self.parent_map[self._node(taxon)]
        out = {'taxon_id': self.taxon_id(taxon)}
        for ancestor in self.ancestors_list(ROOT_RANK, taxon):
            rank = self.rank(ancestor)
            if rank in CANONICAL_RANKS:
                out[rank] = ancestor
        out = {rank: ancestor_rank(rank, taxon)}
        return out

    def phylum(self, taxon, default=None):
        """Return the phylum for the given taxon."""
        return self.ancestor_rank('phylum', taxon, default=default)

    def genus(self, taxon, default=None):
        """Return the genus for the given taxon."""
        return self.ancestor_rank('genus', taxon, default=default)

    def place_microbe(self, taxon):
        """Returns a tuple of (taxonomic id, rank, superkingdom) for the taxon."""
        if taxon == 'root':
            raise TaxonomicRankError('Cannot give superkingdom for root.')
        if self.rank(taxon) in ['subspecies', 'no rank']:
            raise TaxonomicRankError(f'Cannot resolve {taxon} at rank {self.rank(taxon)}.')
        superkingdom = self.ancestor_rank('superkingdom', taxon)
        if superkingdom in ALLOWED_SUPERKINGDOMS:
            return self.taxon_id(taxon), self.rank(taxon), superkingdom
        raise TaxonomicRankError(f'Superkingdom {superkingdom} not allowed.')

    def _tree(self, taxa):
        queue, tree = {self._node(el) for el in taxa}, {}
        root = None
        while len(queue):
            cur_node = queue.pop()
            parent_node = self.parent_map[cur_node]
            if cur_node not in tree:
                tree[cur_node] = {'parent': parent_node, 'children': set()}
            if not parent_node:
                root = cur_node
                continue
            try:
                tree[parent_node]['children'].add(cur_node)
            except KeyError:
                tree[parent_node] = {
                    'parent': self.parent_map[parent_node],
                    'children': set([cur_node])
                }
            queue.add(parent_node)
        return root, tree

    def taxa_sort(self, taxa):
        """Return a list with all elements of taxa in DFS order."""
        taxa, sort = set(taxa), []
        root, tree = self._tree(taxa)

        def dfs(node):
            for child_node in tree[node]['children']:
                child = self._name(child_node)
                if child in taxa:
                    sort.append(child)
                dfs(child_node)
        dfs(root)  # typically 1 is the root node
        return sort

    def all_names(self):
        """Return a list of all scientific names in the tree. Order not guaranteed."""
        return list(self.names_to_nodes.keys())

    @classmethod
    def parse_files(cls, names_filename=None, nodes_filename=None):
        """Return a tree parsed from the given files."""
        names_filename = names_filename if names_filename else environ.get(NAMES_ENV_VAR, NAMES_DEF)
        nodes_filename = nodes_filename if nodes_filename else environ.get(NODES_ENV_VAR, NODES_DEF)
        with gzip.open(names_filename) as names_file, gzip.open(nodes_filename) as nodes_file:
            names_to_nodes, nodes_to_name = {}, {}
            sci_name = list()
            for line in names_file:
                line = line.decode('utf-8')
                tkns = [tkn.strip() for tkn in line.strip().split(NCBI_DELIM)]
                if len(tkns) >= 3:
                    if tkns[3] != 'scientific name':
                        continue
                    node, name = tkns[:2]
                    names_to_nodes[name] = node
                    nodes_to_name[node] = {'name': name, 'rank': None}
                    sci_name.append(name)

            parent_map = {}
            for line in nodes_file:
                line = line.decode('utf-8')
                tkns = [tkn.strip() for tkn in line.strip().split(NCBI_DELIM)]
                if len(tkns) >= 4:
                    node, parent, rank = tkns[:3]
                    if node.isdigit():
                        nodes_to_name[node]['rank'] = rank
                    if node == parent:  # NCBI has a self loop at the root
                        parent = None
                    parent_map[node] = parent
        return cls(parent_map, names_to_nodes, nodes_to_name)
