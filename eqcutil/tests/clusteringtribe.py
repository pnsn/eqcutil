import unittest, os
from pathlib import Path

import numpy as np
import pandas as pd
from obspy.clients.fdsn import Client
from eqcorrscan import Tribe, Template
from eqcorrscan.utils.clustering import space_cluster, space_time_cluster, cluster

from eqcutil.augment.template import rename_templates
from eqcutil.core.clusteringtribe import ClusteringTribe

class TestClusteringTribe(unittest.TestCase):
    cwd = Path(__file__).parent
    TRIBE = Tribe().read(os.path.join(cwd, 'basic_tribe.tgz'))
    TRIBE = rename_templates(TRIBE)

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
        self.assertIn(self.ctribe[0].name, self.ctribe.clusters.index.values)

    def test_add_template_duplicate(self):
        self.ctribe.add_template(self.tribe[0].copy())
        with self.assertRaises(AttributeError):
            self.ctribe.add_template(self.tribe[0].copy(), rename_duplicates=False)
        self.ctribe.add_template(self.tribe[0], rename_duplicates=True)
        self.assertEqual(len(self.ctribe), 2)
        self.assertEqual(self.ctribe[0].name + '__0', self.ctribe[1].name)
    
    def test_extend_tribe(self):
        self.ctribe.extend(self.tribe)
        self.assertEqual(len(self.ctribe), len(self.tribe))
        self.assertEqual(self.ctribe.templates, self.tribe.templates)
        self.assertTrue(all(_t.name in self.ctribe.clusters.index.values) for _t in self.ctribe.templates)

    def test_space_cluster(self):
        self.ctribe.extend(self.tribe.copy())
        kwargs = {'d_thresh': 3., 'show': False}
        # Assert just id_no column in clusters
        self.assertEqual(set(self.ctribe.clusters.columns.values), set(['id_no'])) 
        # Run original method
        tribes = self.tribe.copy().cluster('space_cluster', **kwargs)
        ## TEST SPACE_CLUSTER
        self.ctribe.cluster('space_cluster', **kwargs)
        # Assert increase in clusters size
        self.assertEqual(self.ctribe.clusters.shape, (len(self.ctribe), 2))
        # Test that new implementation matches function-oriented original
        for name, row in self.ctribe.clusters.iterrows():
            idx = row.space_cluster
            self.assertIn(self.ctribe.select(name), tribes[idx].templates)
        self.assertEqual(self.ctribe.cluster_kwargs['space_cluster'], kwargs)

    def test_space_time_cluster(self):
        self.ctribe.extend(self.tribe.copy())
        kwargs = {'d_thresh': 5., 't_thresh': 86400}
        # Assert just id_no column in clusters
        self.assertEqual(set(self.ctribe.clusters.columns.values), set(['id_no'])) 
        # Run original method
        tribes = self.tribe.copy().cluster('space_time_cluster', **kwargs)
        # Run new implementation
        self.ctribe.cluster('space_time_cluster', **kwargs)
        # Assert that clusters has increased by one column
        self.assertEqual(self.ctribe.clusters.shape, (len(self.ctribe), 2))
        # Assert that the column names are correct
        self.assertEqual(set(self.ctribe.clusters.columns.values), set(['id_no','space_time_cluster']))
        # Check that results are the same
        for name, row in self.ctribe.clusters.iterrows():
            idx = row.space_time_cluster
            self.assertIn(self.ctribe.select(name), tribes[idx].templates)
        # Check that cluster_kwargs matches the kwargs
        self.assertEqual(self.ctribe.cluster_kwargs['space_time_cluster'], kwargs)

    def test_correlation_cluster(self):
        self.ctribe.extend(self.tribe.copy())
        self.assertIsNone(self.ctribe.dist_mat)
        kwargs = {'show': False, 'corr_thresh': 0.3, 'allow_individual_trace_shifts': False,
                  'save_corrmat': True, 'replace_nan_distances_with': 1, 'cores': 'all'}
        # Run initial implementation 
        groups = cluster([(_t.st, _e) for _e, _t in enumerate(self.tribe.copy())], **kwargs)
        # Run new implementation
        self.ctribe.cluster('correlation_cluster', **kwargs)
        # Assert columns in ctribe.clusters
        self.assertEqual(set(self.ctribe.clusters.columns.values), set(['id_no','correlation_cluster']))
        # Assert shape of clusters dataframe
        self.assertEqual(self.ctribe.clusters.shape, (len(self.ctribe), 2))
        # Assert that results are the same
        for name, row in self.ctribe.clusters.iterrows():
            id_no = row.id_no
            grp_no = row.correlation_cluster
            self.assertIn(id_no, [entry[1] for entry in groups[grp_no]])
        # Assert that dist_mat populated
        self.assertIsInstance(self.ctribe.dist_mat, np.ndarray)
        # Assert that dist_mat is scaled correctly
        self.assertEqual(self.ctribe.dist_mat.shape, (len(self.ctribe), len(self.ctribe)))
        # Assert that kwargs are saved
        self.assertEqual(self.ctribe.cluster_kwargs['correlation_cluster'], kwargs)
        os.remove(Path().cwd()/'dist_mat.npy')



