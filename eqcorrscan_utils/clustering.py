import os
from pathlib import Path

from eqcorrscan import Template, Tribe

import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform
from eqcorrscan.utils.clustering import handle_distmat_nans


def catalog_cluster(catalog, thresh, metric='distance', show=False):
    """Alias to EQcorrscan :meth:`~eqcorrscan.utils.clustering.catalog_cluster`
    that does spatial or clustering on events in an ObsPy :class:`~obspy.core.event.Catalog`
    object. 

    https://eqcorrscan.readthedocs.io/en/latest/submodules/autogen/eqcorrscan.utils.clustering.catalog_cluster.html#eqcorrscan.utils.clustering.catalog_cluster

    :param catalog: catalog with events to cluster
    :type catalog: :class:`~obspy.core.event.Catalog`
    :param thresh: threshold value (km for *metric="distance"*, sec for *metric="time"*)
    :type thresh: float
    :param metric: threshold metric, either 'distance' or 'time', defaults to 'distance'
    :type metric: str, optional
    :param show: should the clusters be plotted? Defaults to False
    :type show: bool, optional
    :return:
     - **clusters** (*list* of *obspy.core.event.Catalog*) - sub-catalogs containing
        clustered events.
    """    
    clusters = catalog_cluster(catalog, thresh, metric=metric, show=show)
    return clusters


def tribe_cluster(tribe, save_corrmat=os.path.join('.','tribe_clustering_outputs'), **kwargs):
    if not isinstance(tribe, Tribe):
        raise TypeError(f'tribe must be type eqcorrscan.core.match_filter.tribe.Tribe')
    
    template_list = _compose_template_list(tribe)
    # Force save
    if save_corrmat:
        savename = kwargs['save_corrmat'].copy()
        kwargs.update({'save_corrmat': True})

    groups = cluster(template_list, **kwargs)
    



def _compose_template_list(tribe):
    """Compose a template list formatted as an input
    for :meth:`~eqcorrscan.util.clustering.cluster`

    template_list has the format:
     [(template.st, template.name), ...]
     where stream is template.st, and index is a unique
     integer key.

     i.e., a list of 2-tuples each containing: 
            [0] the template waveforms (`template.st`)
            [1] the template name (`template.name`)
    
    :param tribe: tribe containing templates
    :type tribe: eqcorrscan.core.match_filter.tribe.Tribe
    :return: tlist
    :rtype: list
    """    
    tlist = [(template.st, template.name) for template in tribe]
    return tlist


def _save_template_clustering_output(save_dir, groups, savename='clusters'):

    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    with open(os.path.join(save_dir, f'{savename}.csv'), 'w') as ofile:
        ofile.write('template_name,group_id\n')
        for _e, group in enumerate(groups):
            for entry in group:
                ofile.write(f'{entry[1]},{_e}\n')

    src = os.path.join(Path().cwd(), 'dist_mat.npy')
    dest = os.path.join(save_dir, f'{savename}_dist_mat.npy')
    if os.path.isfile(src):
        os.rename(src,dest)

def reconstitute_linkage(distance_matrix_file, fill_value=1, **kwargs):
    """Re-load a saved `dist_mat.npy`-type file saved by
    :meth:`~eqcorrscan.util.cluster.clustering` and re-create
    the `Z` input to :meth:`~scipy.cluster.dendrogram` via the
    :meth:`~scipy.cluster.hierarchical.linkage` method and provide
    an option for a NaN fill_value

    **fill_value** corresponds to the `replace_nan_distances_with` parameter
    in :meth:`~eqcorrscan.util.cluster.clustering`.
    
    .. rubric:: Notes on `fill_value`
        Useful values for fill_value might include:
        * 1 - treat all missing channel-wise correlations as
            contributing the maximum distance (i.e., 0 correlation)
            to the template-averaged correlation / distance.
            This produces a conservative estimate of template
            correlation strength as it down-weights the importance
            of event pairs that have missing data.
        * 'mean' - fill missing values with the average
            correlation / distance value calculated from
            non-NaN elements. This results in not
            down-weighting a correlation due to missing data.

    :param distance_matrix_file: name of the numpy distance matrix
        file to load
    :type distance_matrix_file: str
    :param corr_thresh: correlation coefficient threshold for
        segmenting clusters, defaults to 0.4
    :type corr_thresh: float, optional
    :param fill_value: fill value for NaN elements of the
        distance matrix, defaults to 1.
        see note above
    :type fill_value: int, optional
    :return: linkage_array - output o
    :rtype: scipy.cluster.hierarchy.linkage
    """    
    dist_mat = np.load(distance_matrix_file)
    dist_mat = handle_distmat_nans(dist_mat, replace_nan_distances_with=fill_value)
    dist_vect = squareform(dist_mat)
    Z = linkage(dist_vect, **kwargs)

    return Z

def make_dendrogram(dist_file, corr_threshold=0.4, fill_value=1, linkage_kwargs={}, **kwargs):
    """Re-create a dendrogram figure from the dist_mat.npy
    file saved by :meth:`~eqcorrscan.utils.cluster.clustering`

    :param dist_file: distance matrix file to load
    :type dist_file: str
    :param corr_threshold: correlation coefficient threshold,
        for delineating clusters, defaults to 0.4
    :type corr_threshold: float, optional
    :param fill_value: NaN fill value for the distance matrix,
         defaults to 1. 
         also see :meth:`~eqcorrscan_utils.methods.clustering.reconstitute_linkage`
    :type fill_value: int, optional
    :param linkage_kwargs: Additional key-word arguments to pass,
        to `reconstitute_linkage`, defaults to {}
    :type linkage_kwargs: dict, optional
    :param kwargs: key-word argument collector to pass to
        :meth:`~scipy.cluster.dendrogram`

    :return: dictionary of structures computed to render the dendrogram.
        default output of :meth:`~scipy.cluster.dendrogram`
    :rtype: dict
    """    
    Z = reconstitute_linkage(dist_file, fill_value=fill_value, **linkage_kwargs)

    R = dendrogram(Z, color_threshold=1 - corr_threshold,
               distance_sort='ascending', ax=None)
    return R