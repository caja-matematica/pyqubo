import unittest

from pyqubo.utils.ising_formulations import *
from pyqubo.utils.graph_converter import *

import networkx as nx


class TestGraphConverter(unittest.TestCase):

    def test_graph_to_qbits(self):
        G = nx.Graph()
        G.add_edges_from(
            [(0, 1), (0, 2), (1, 3), (2, 3), (3, 5), (4, 5), (3, 6), (4, 7), (5, 6), (5, 7), (6, 8), (7, 8), (8, 9)])

        G_copy, node_hash = graph_to_qbits(G)

        self.assertEqual(len(node_hash), len(G_copy))
        self.assertTrue('label' in G_copy.node[0])
        self.assertEqual(G_copy.node[0]['label'], 'x_0')


class TestIsingFormulations(unittest.TestCase):

    def test_vertex_cover(self):
        # Set up a Networkx Graph for a small known problem
        G = nx.Graph()
        G.add_edges_from(
            [(0, 1), (0, 2), (1, 3), (2, 3), (3, 5), (4, 5), (3, 6), (4, 7), (5, 6), (5, 7), (6, 8), (7, 8), (8, 9)])

        q1 = generate_vertex_cover_qubo(G)

        known_qubo = {(3, 3): -35.0,
                     (5, 5): -35.0,
                     (6, 6): -25.0,
                     (7, 7): -25.0,
                     (8, 8): -25.0,
                     (0, 0): -15.0,
                     (1, 1): -15.0,
                     (2, 2): -15.0,
                     (4, 4): -15.0,
                     (9, 9): -5.0,
                     (0, 1): 10.0,
                     (0, 2): 10.0,
                     (1, 3): 10.0,
                     (2, 3): 10.0,
                     (3, 5): 10.0,
                     (3, 6): 10.0,
                     (4, 5): 10.0,
                     (4, 7): 10.0,
                     (5, 6): 10.0,
                     (5, 7): 10.0,
                     (6, 8): 10.0,
                     (7, 8): 10.0,
                     (8, 9): 10.0}
        const = 130.0
        converted_to_pyq = {("x_{}".format(u), "x_{}".format(v)): known_qubo[u,v] for u,v in known_qubo}
        pyq, pyqconst = q1.compile().to_qubo()

        self.assertDictEqual(converted_to_pyq, pyq)
        self.assertEqual(const, pyqconst)


if __name__ == "__main__":
    unittest.main()