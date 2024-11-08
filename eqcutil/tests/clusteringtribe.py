import unittest, os
from pathlib import Path

import pandas as pd
from obspy.clients.fdsn import Client
from eqcorrscan import Tribe, Template

from eqcutil.core.clusteringtribe import ClusteringTribe

class TestClusteringTribe(unittest.TestCase):
    cwd = Path(__file__).parent
    TRIBE = Tribe().read(os.path.join(cwd, 'basic_tribe.tgz'))

    def setUp(self):
        self.ctribe = ClusteringTribe()
        self.tribe = self.TRIBE.copy()
    
    def tearDown(self):
        del self.ctribe
        del self.tribe


    def test___init__(self):
        self.assertIsInstance(self.ctribe, Tribe)
        self.assertIsInstance(self.ctribe, ClusteringTribe)
        self.assertIsInstance(self.ctribe.clusters, pd.DataFrame)
        self.assertIsNone(self.ctribe.dist_mat)
        self.assertIsInstance(self.ctribe.cluster_kwargs, dict)

    def test_add_template(self):
        self.ctribe.add_template(self.tribe.templates[0])
        self.assertEqual(len(self.ctribe), 1)
        self.assertIsInstance(self.ctribe[0], Template)
        self.assertIn(self.ctribe[0].name, self.ctribe.clusters.name.values)

    def test_add_template_duplicate(self):
        self.ctribe.add_template(self.tribe[0])
        with self.assertRaises(AttributeError):
            self.ctribe.add_template(self.tribe[0], rename_duplicates=False)