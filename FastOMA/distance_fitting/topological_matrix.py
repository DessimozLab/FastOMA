"""
OMA HOG distance fitting
Fits distances from a subset of trees to OMA HOGs.

(C) 2025 Comparative Genomics Group <contact@omabrowser.org>

This file is part of hog_distances.

hog_distances is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

hog_distances is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with hog_distances.  If not, see <http://www.gnu.org/licenses/>.
"""

from datetime import datetime
from scipy.sparse import csr_matrix
import numpy as np

from .nngs import nngs


class TopologicalMatrix(object):
    def __init__(self, t):
        self.t = t
        self.t.encode_bipartitions(collapse_unrooted_basal_bifurcation=False)
        self.n = len(self.t.leaf_nodes())
        self._set_edges()
        
        t_leaf_taxa = set(map(lambda n: n.taxon, t.leaf_node_iter()))
        self._leaf_edge_tab = list(map(lambda x: self._get_taxa_edges(x.label),
                                       filter(lambda x: x in t_leaf_taxa,
                                              self.t.taxon_namespace)))

        conv_to_int = lambda x: int(x) if x.isdigit() else x
        self.taxa_map = dict(map(lambda x: (conv_to_int(x[1].label), x[0]),
                                 enumerate(filter(lambda x: x in t_leaf_taxa,
                                                  self.t.taxon_namespace))))
        #self.taxa = list(map(int, self.t.taxon_namespace.labels()))

    def _get_taxa_edges(self, x):
        edges = set()
        n = self.t.find_node_with_taxon_label(x)
        while n is not self.t.seed_node:
            if n.edge.__i is not None:
                edges.add(n.edge.__i)
            n = n.parent_node
        return edges

    def _set_edges(self):
        self.edges = self.t.leaf_edges()
        self.edges += list(self.t.levelorder_edge_iter(lambda e:
                                                       e.is_internal()))[::-1]

        self.n_edges = len(self.edges) - 1
        for (i, e) in enumerate(self.edges):
            e.__i = i
        if len(self.t.seed_node.child_nodes()) == 2:
            self.t.seed_node.child_edges()[0].__i = min(map(lambda e: e.__i, self.t._seed_node.child_edges()))
            self.t.seed_node.child_edges()[1].__i = None
            #self.edges[-2].__i = None
            self.n_edges -= 1
        self.edges[-1].__i = None

    def get_edges(self, i, j):
        return np.array(sorted(self._leaf_edge_tab[i] ^ self._leaf_edge_tab[j]),
                        dtype=np.uint32)

    def format(self, pairs):
        n = len(pairs)
        #fixed_d = []
        #edges_fixed = set()
        #fixed_edge_ii = []

        #indptr = [0]  #np.zeros(n+1, dtype=np.int64)
        indptr = np.zeros(n+1, dtype=np.int64)
        indices = []
        ind_count = 0

        max_row = 0

        ## include distances on the tree already.
        #for e in self.t.edges():
        #    if e.__i is not None and e.length is not None:
        #        if e.is_internal():
        #            continue
        #        edges_fixed.add(e.__i)
        #        fixed_edge_ii.append(e.__i)
        #        fixed_d.append(e.length)
        #        indices.append(np.array([e.__i]))
        #        ind_count += 1
        #        indptr.append(ind_count)

        #pair_mask = [False] * len(pairs)
        for (i, p) in enumerate(pairs):
            z = self.get_edges(*map(self.taxa_map.__getitem__, p))
            indices.append(z)
            ind_count += len(z)
            indptr[i+1] = ind_count

            max_row = max(max_row, len(z))

            #if len(set(z) - edges_fixed) > 0:
            #    # observation over edge without fixed distance
            #    pair_mask[i] = True
            #    indices.append(z)
            #    ind_count += len(z)
            #    indptr.append(ind_count)

        # form M
        indices = np.concatenate(indices)
        data = np.ones(len(indices), dtype=np.int64)  # we have overflow potential with scalars in the numba code
        return csr_matrix((data, indices, indptr)).tocsc()  # More efficient for computing A.T@A
        #return (csr_matrix((data, indices, indptr)).tocsc(),  # More efficient for computing A.T@A
        #        np.array(fixed_d),
        #        pair_mask,
        #        np.array(fixed_edge_ii))

    def solve(self, pairs, d): #, fit_missing_free):  #, convergence_type=2, initial=1):
        if self.n > 2:
            #(M, fixed_d, pair_mask, edges_fixed) = self.format(pairs)
            M = self.format(pairs)
            #d1 = np.concatenate((fixed_d, d[pair_mask]))

            t0 = datetime.now()
            #(x, it, r) = nngs(M, d, converge_type=convergence_type, initial=initial)
            #(x, it, r) = nngs(M, d1, edges_fixed, fit_missing_free)
            (x, it, r) = nngs(M, d)
            t = (datetime.now() - t0).total_seconds()
        else:
            (x, it, r) = (np.array([np.mean(d),]), 0, 0)
            t = 0.0

        # fit distances back onto tree
        self.fit_distances(x)
        return (it, r, t)

    def fit_distances(self, x):
        for e in filter(lambda e: (e.__i is not None),
                        self.t.edges()):
            e.length = x[e.__i]

        if len(self.t.seed_node.child_edges()) == 2:
            # allow for bifurcation at root of hog tree
            hog_root_e1 = self.t.seed_node.child_edges()[0]
            hog_root_e2 = self.t.seed_node.child_edges()[1]
            if hog_root_e1.length is not None:
                hog_root_e1.length = hog_root_e2.length = (hog_root_e1.length / 2)

        return self.t
