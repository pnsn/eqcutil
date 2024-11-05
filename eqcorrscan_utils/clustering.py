"""
:module: eqc_utils.clustering
:compiler: Nathan T. Stevens
:email: ntsteven@uw.edu
:org: Pacific Northwest Seismic Network
:license: GNU GPLv3
:purpose:
    This module provides a "one-stop-shop" for assorted clustering
    methods from EQcorrscan and helper methods facilitating their
    use.
:attribution:
    Any use of these methods should include an ackowledement of
    the EQcorrscan project or paper (e.g., Chamberlain et al., 2017)
"""
import os, logging
from pathlib import Path
from multiprocessing import cpu_count

import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform

from eqcorrscan import Template, Tribe
from eqcorrscan.utils.clustering import handle_distmat_nans, distance_matrix, fcluster



Logger = logging.getLogger(__name__)

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


def cross_corr_cluster(
        templates, show=False, corr_thresh=0.3, shift_len=0,
        allow_individual_trace_shifts=True, fill_value=None,
        cores='all', save_path=False, save_subtribes=True,
        **kwargs):
    """
    Cluster template waveforms based on average correlations.

    Adapted from :meth:`~eqcorrscan.util.clustering.cluster` to
    provide additional file saving/formatting support to make
    it easier to re-analyze tribe clustering analysis outputs.

    Parameters
    ----------

    :type templates: iterable collection of
        :class:`~eqcorrscan.core.match_filter.template.Template` objects
    :param templates:
        list-like
    :type show: bool
    :param show: plot linkage on screen if True, defaults to True
    :type corr_thresh: float
    :param corr_thresh: Cross-channel correlation threshold for grouping
    :type shift_len: float
    :param shift_len: How many seconds to allow the templates to shift
    :type allow_individual_trace_shifts: bool
    :param allow_individual_trace_shifts:
        Controls whether templates are shifted by shift_len in relation to the
        picks as a whole, or whether each trace can be shifted individually.
        Defaults to True.
    :type save_path: bool or path-like str
    :param save_path: If False, will not save, otherwise the
        distance matrix, shift matrix, input parameters, and template
        name index are saved to the specified directory. Directory is
        automatically generated if it doesn not exist
    :type fill_value: None, 'mean', 'min', or float
    :param fill_value:
        Controls how the clustering handles nan-distances in the distance
        matrix. None/False only performs a check, while other choices (e.g.,
        1, 'mean', 'min' or float) replace nans in the distance matrix.
    :type cores: int
    :param cores:
        number of cores to use when computing the distance matrix, defaults to
        'all' which will work out how many cpus are available and hog them.

    :returns:
        List of groups. Each group is a list of
        :class:`obspy.core.stream.Stream` making up that group.

    Original Documentation from `cluster`
    ------------------------------------

    Function to take a set of templates and cluster them, will return groups
    as lists of streams.  Clustering is done by computing the cross-channel
    correlation sum of each stream in stream_list with every other stream in
    the list.  :mod:`scipy.cluster.hierarchy` functions are then used to
    compute the complete distance matrix, where distance is 1 minus the
    normalised cross-correlation sum such that larger distances are less
    similar events.  Groups are then created by clustering the distance matrix
    at distances less than 1 - corr_thresh.

    When distance_matrix contains NaNs (event pairs that cannot be directly
    compared), then the mean correlation between templates is used instead of
    NaN (see https://github.com/eqcorrscan/EQcorrscan/issues/484).

    Will compute the distance matrix in parallel, using all available cores.
    The method, metric, and order to compute linkage from the distance matrix
    can be controled with parameters from scipy.cluster.hierarchy.linkage as
    kwargs.

    Modifications
    -------------
    
    **templates** -- Supercedes **template_list**. Provides the expected
    "list of tuples" format for sub-methods, creating tuples as
    (template.st, template.name).

    **save_path** -- This method now saves the distance matrix, the shift
    matrix, and the shift dictionary output by sub-method
    :meth:`~eqcorrscan.util.clustering.distance_matrix` to a specified
    save directory. Supercedes the **save_corrmat** argument.

    **fill_value** -- shorthand alias for "replace_nan_distances_with"

    
    """
    if cores == 'all':
        num_cores = cpu_count()
    else:
        num_cores = cores

    if isinstance(save_path, str):
        if not os.path.exists(save_path):
            os.makedirs(save_path)
        # Save parameters
        with open(os.path.join(save_path, 'params.csv'), 'w') as _f:
            _f.write(f'show,{show}\n')
            _f.write(f'corr_thresh,{corr_thresh}\n')
            _f.write(f'shift_len,{shift_len}\n')
            _f.write(f'allow_individual_trace_shifts,{allow_individual_trace_shifts}\n')
            _f.write(f'save_path,{save_path}\n')
            _f.write(f'fill_value,{fill_value}\n')
            _f.write(f'cores,{cores}\n')
            for _k, _v in kwargs.items():
                _f.write(f'{_k},{_v}\n')
    
    if 'save_corrmat' in kwargs.keys():
        raise DeprecationWarning('save_corrmat has been replaced with save_path for this version of `correlate`')

    template_list = _compose_template_list(templates)

    # Save the index/name relationship
    if save_path:
        with open(os.path.join(save_path,'index.csv'),'w') as _f:
            _f.write(f'index, name\n')
            for _e, temp in enumerate(templates):
                _f.write(f'{_e},{temp.name}\n')
        
    # Extract only the Streams from stream_list
    stream_list = [x[0] for x in template_list]

    # Compute the distance matrix
    Logger.info('Computing the distance matrix using %i cores' % num_cores)
    dist_mat, shift_mat, shift_dict = distance_matrix(
        stream_list=stream_list,
        shift_len=shift_len,
        cores=num_cores,
        replace_nan_distances_with=fill_value,
        allow_individual_trace_shifts=allow_individual_trace_shifts)
    # Save outputs
    if save_path:
        # Save distance matrix
        np.save(os.path.join(save_path, 'dist_mat.npy'), dist_mat)
        Logger.info('Saved the distance matrix as dist_mat.npy')
        np.save(os.path.join(save_path, 'shift_mat.npy') shift_mat)
        Logger.info('Saved the shift matrix as shift_mat.npy')
    
    # Calculate linkage
    Logger.info('Computing linkage')
    dist_mat = handle_distmat_nans(dist_mat,
                                   replace_nan_distances_with=fill_value)
    dist_vec = squareform(dist_mat)
    Z = linkage(dist_vec, **kwargs)

    # Get the indices of the groups
    Logger.info('Clustering')
    indices = fcluster(Z, t=1 - corr_thresh, criterion='distance')
    # Indices start at 1...
    group_ids = list(set(indices))  # Unique list of group ids
    Logger.info(' '.join(['Found', str(len(group_ids)), 'groups']))
    # Convert to tuple of (group id, stream id)
    indices = [(indices[i], i) for i in range(len(indices))]

    # Sort by group id
    indices.sort(key=lambda tup: tup[0])
    groups = []
    Logger.info('Extracting and grouping')
    for group_id in group_ids:
        group = []
        for ind in indices:
            if ind[0] == group_id:
                group.append(template_list[ind[1]])
            elif ind[0] > group_id:
                # Because we have sorted by group id, when the index is greater
                # than the group_id we can break the inner loop.
                # Patch applied by CJC 05/11/2015
                groups.append(group)
                break
    # Catch the final group
    groups.append(group)

    # OPTIONAL - Save the groups as sub-tribes
    if save_subtribes and save_path:
        isave = os.path.join(save_path,'sub_tribes')
        if not os.path.exists(isave):
            os.makedirs(isave)
        for _e, group in enumerate(groups):
            itribe = Tribe(templates=group)
            itribe.write(os.path.join(isave,f'sub_tribe_{_e}.tgz'))
            isave = os.path.join(save_path, f'group_{_e}')

    # OPTIONAL - Show the dendrogram
    if show:
        Logger.info('Plotting the dendrogram')
        dendrogram(Z, color_threshold=1 - corr_thresh,
                   distance_sort='ascending')
        plt.show()
    return groups


def tribe_cluster(templates, save_path='./cluster_results',
                  corr_thresh=0.3, shift_len=1., fill_value = 1, 
                  allow_individual_trace_shifts=False, cores='all', show=True, method='linear',
                  metrci='euclidian', optimal_ordering=False):
    """Conduct waveform-based clustering of :class:`~eqcorrscan.core.match_filter.template.Template` objects
    using the :meth:`~eqcorrscan.util.clustering.cluster` method and record the parameter inputs
    and results (including the distance matrix) in a specified directory.

    We provide an extended set of input parameters to encompass patameters for this method,
    and the underlying :meth:`~eqcorrscan.util.clustering.cluster` and :meth:`~scipy.cluster.hierarchy.linkage`
    methods. 

    `tribe_cluster` Parameters
    --------------------------
    :param templates: iterable collection of template objects
    :type: eqcorrscan.core.match_filter.tribe.Tribe or similar
    :param save_path: relative or absolute path to the directory in which
        clustering results should be saved, defaults to './clustering_results'.
        This directory path will be generated if it does not already exist
    :type save_path: str, optional

    
    `cluster` Parameters
    --------------------

    :param corr_thresh: cross-correlation threshold for delineating
        template clusters, defaults to 0.3
    :type corr_thresh: float, optional
    :param shift_len: maximum template shift in seconds to seek correlation maxima,
        defaults to 1.
    :type shift_len: float, optional
    :param fill_value: alias for "replace_nan_distances_with" parameter in
        :meth:`~eqcorrscan.util.clustering.cluster`, defaults to 1.
        See notes below
    :type fill_value: scalar, str, or NoneType, optional
    :param allow_individual_trace_shifts: Controls wheterh templates are shifted
        by **shift_len** in relation to the picks as a whole, or
        whether they can be shifted individually. Defaults to False

    `linkage` Parameters
    --------------------
    :param method: method to use for linkage


    Files Generated
    ---------------
     - **params.csv** -- all input parameters to this method as a CSV
        with parameter names in the left column and values in the right column
     - **groups_cct###.csv** -- names of templates (left column) and their group number (right column)
        based on the provided **corr_thresh** value (_cct### = _cct{corr_thresh: %.2f})
     - **dist_mat.npy*** -- distance matrix saved by 

    """    
    # Must be an iterable composed only of Template objects
    if not hasattr(templates, '__iter__'):
        raise AttributeError(f'input "templates" is not iterable')
    elif not all(isinstance(_e, Template) for _e in templates):
        raise TypeError(f'not all elements of iterable "templates" are type eqcorrscan.core.match_filter.template.Template')
    else:
        pass
    
    if not isinstance(save_path, str):
        raise TypeError
    elif not os.path.exists(save_path):
        os.makedirs(save_path)
    else:
        pass

    if fill_value not in ['mean','max']

    template_list = _compose_template_list(templates)

    # Force save
    if save_corrmat:
        savename = kwargs['save_corrmat'].copy()
        kwargs.update({'save_corrmat': True})

    groups = cluster(template_list, **kwargs)
    

def _compose_template_list(templates):
    """Compose a template list formatted as an input
    for :meth:`~eqcorrscan.util.clustering.cluster`

    template_list has the format:
     [(template.st, template.name), ...]
     where stream is template.st, and index is a unique
     integer key.

     i.e., a list of 2-tuples each containing: 
            [0] the template waveforms (`template.st`)
            [1] the template name (`template.name`)
    
    :param templates: iterable collection of template objects
    :type tribe: iterable that yields eqcorrscan.core.match_filter.tribe.Tribe elements
    :return: 
     - **tlist** (*list*) -- list of tuples with (stream, name) structure
    """
    if not hasattr(templates, '__iter__'):
        raise AttributeError('"templates" is not iterable')
    elif not all(isinstance(_e, Template) for _e in templates):
        raise TypeError('not all elements of "templates" are type eqcorrscan.core.match_filter.template.Template')
    else:
        pass
    tlist = [(template.st, template.name) for template in templates]
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
