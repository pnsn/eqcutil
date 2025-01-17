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

TODO: Merge cross_corr_cluster and tribe_cluster
"""
import os, logging
from pathlib import Path
from multiprocessing import cpu_count

from tqdm import tqdm
import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform

from obspy import Trace, Stream
from eqcorrscan import Template, Tribe
from eqcorrscan.utils.clustering import cluster, handle_distmat_nans, cross_chan_correlation

from eqcutil.util.logging import rich_error_message


module_logger = logging.getLogger(__name__)


def distance_matrix(
        stream_list, shift_len=0.0,
        replace_nan_distances_with='mean',
        allow_individual_trace_shifts=True,
        cores=1, progress_bar=True
    ):
    """ :meth:`~eqcorrscan.utils.clustering.distance_matrix` with 
    a few additional safety catches added, inclusion of a tqdm progress
    bar, and modified default input values

Compute distance matrix for waveforms based on cross-correlations.

    Function to compute the distance matrix for all templates - will give
    distance as 1-abs(cccoh), e.g. a well correlated pair of templates will
    have small distances, and an equally well correlated reverse image will
    have the same distance as a positively correlated image - this is an issue.

    :type stream_list: list
    :param stream_list:
        List of the :class:`obspy.core.stream.Stream` to compute the distance
        matrix for
    :type shift_len: float, optional
    :param shift_len: How many seconds for templates to shift
        Defaults to 0.
    :type allow_individual_trace_shifts: bool, optional
    :param allow_individual_trace_shifts:
        Controls whether templates are shifted by shift_len in relation to the
        picks as a whole, or whether each trace can be shifted individually.
        Defaults to True.
    :type replace_nan_distances_with: None, 'mean', 'min', or float, optional
    :param replace_nan_distances_with:
        Controls how the clustering handles nan-distances in the distance
        matrix. None/False only performs a check, while other choices (e.g.,
        1, 'mean', 'min' or float) replace nans in the distance matrix.
        Defaults to 'mean'
    :type cores: int
    :param cores: Number of cores to parallel process using, defaults to 1.

    :returns:
        - distance matrix (:py:class:`numpy.ndarray`) of size
          len(stream_list)**2
        - shift matrix (:py:class:`numpy.ndarray`) containing shifts between
          traces of the sorted streams. Size is len(stream_list)**2 * x, where
          x is 1 for shift_len=0 and/or allow_individual_trace_shifts=False.
          Missing correlations are indicated by nans.
        - shift dict (:py:class:`dict`):
          dictionary of (template_id: trace_dict) where trace_dict contains
          (trace.id: shift matrix (size `len(stream_list)**2`) for trace.id)

    .. warning::
        Because distance is given as :math:`1-abs(coherence)`, negatively
        correlated and positively correlated objects are given the same
        distance.

    .. note::
        Requires all traces to have the same sampling rate and same length.
    """    
    n_streams = len(stream_list)
    # May have to allow duplicate channels for P- and S-picks at each station
    # stream_list = [st.sort() for st in stream_list]
    # uniq_traces = set([tr.id for st in stream_list for tr in st])
    # n_uniq_traces = len(uniq_traces)

    # Added safety catch for stream_list checking if there are multiple traces of the same id
    uniq_traces = set([])
    for st in stream_list:
        if not isinstance(st, Stream):
            module_logger.critical('One or more input elements of stream_list are not type obspy.Stream')
        else:
            st.sort()
        utid = set([tr.id for tr in st])
        if len(utid) != len(st):
            module_logger.warning('Multiple traces with same ID - may result in an abberent behavior when populating the shift-matrix')
            module_logger.warning('Consider using stream.merge() on input traces before using this method.')
        uniq_traces = uniq_traces.union(utid)

    n_uniq_traces = len(uniq_traces)

    # Initialize square matrix
    dist_mat = np.zeros([n_streams, n_streams])
    shift_mat = np.empty([n_streams, n_streams, n_uniq_traces])
    shift_mat[:] = np.nan
    shift_dict = dict()
    i = -1
    if len(stream_list) == 0:
        breakpoint()
    for master in tqdm(stream_list, disable=not progress_bar):
        i += 1
        #module_logger.debug(f'Distance matrix with master {i+1} of {len(stream_list)}')
        dist_list, shift_list = cross_chan_correlation(
            st1=master, streams=stream_list, shift_len=shift_len,
            allow_individual_trace_shifts=allow_individual_trace_shifts,
            xcorr_func='fftw', cores=cores)
        dist_mat[i] = 1 - dist_list
        master_ids = [tr.id for tr in master]
        master_trace_indcs = [
            j for j, tr_id in enumerate(uniq_traces) if tr_id in master_ids]
        # Sort computed shifts into shift-matrix. shift_list could contain a
        # nan-column that needs to be ignored here (only when earliest trace is
        # missing)
        # Wrapped this segment to avoid crashes coming from multiple traces of the same id
        try:
            shift_mat[np.ix_([i], list(range(n_streams)), master_trace_indcs)] = (
                shift_list[:, ~np.all(np.isnan(shift_list), axis=0)])
        except Exception as e:
            module_logger.warning(rich_error_message(e))
            module_logger.warning(f'Abberent behavior arose while populating shift-matrix for tempate index {i} - leaving as NaN entry')

            
        # Add trace-id with corresponding shift-matrix to shift-dictionary
        shift_mat_list = [shift_mat[:, :, mti] for mti in master_trace_indcs]
        trace_shift_dict = dict(zip(master_ids, shift_mat_list))
        shift_dict[i] = trace_shift_dict
    if shift_len == 0:
        dist_mat = handle_distmat_nans(
            dist_mat, replace_nan_distances_with=replace_nan_distances_with)
    else:
        # get the shortest distance for each correlation pair
        dist_mat_shortest = np.minimum(dist_mat, dist_mat.T)
        # Indicator says which matrix has shortest dist: value 0: mat2; 1: mat1
        mat_indicator = dist_mat_shortest == dist_mat
        mat_indicator = np.repeat(mat_indicator[:, :, np.newaxis],
                                  n_uniq_traces, axis=2)[:, :]
        # Get shift for the shortest distances
        shift_mat = (
            shift_mat * mat_indicator +
            np.transpose(shift_mat, [1, 0, 2]) * (1 - mat_indicator))
        dist_mat = dist_mat_shortest
    # Squeeze matrix to 2 axis (ignore nans) if 3rd dimension not needed
    if shift_len == 0 or allow_individual_trace_shifts is False:
        shift_mat = np.nanmean(shift_mat, axis=2)
    np.fill_diagonal(dist_mat, 0)
    return dist_mat, shift_mat.squeeze(), shift_dict


def compute_pariwise_cross_correlations(
        template_list,
        shift_len=0,
        allow_individual_trace_shifts=True,
        replace_nan_distances_with='mean',
        cores='all'):
    """
    The first half of :meth:`~eqcorrscan.utils.clustering.cluster`
    allowing exposure of the dist_mat and shift_mat for subsequent
    (re)analyses without having to re-run the cross correlation analysis

    :param template_list: 
    
    """
    if cores == 'all':
        num_cores = cpu_count()
    else:
        num_cores = cores
    # Extract streams from stream/index tuples
    stream_list = []
    for _x in template_list:
        if len(_x) != 2:
            module_logger.critical('template list entries must be 2-tuples of type (obspy.Stream, int)')
        if not isinstance(_x[1], int):
            module_logger.warning('template list index is not type int, may result in abberent behavior')
        # Allow some resilliance if a trace is passed instead of a stream
        if not isinstance(_x[0], Stream):
            if isinstance(_x[0], Trace):
                module_logger.debug(f'Trace for "{_x[0].id}" passed as template - wrapping in an obspy.Stream')
                stream_list.append(Stream(_x[0]))
            else:
                module_logger.critical('One or more elements in template_list is not type obspy.Stream')
        else:
            stream_list.append(_x[0])
    # Compute distance matrix, shift matrix, and shift dictionary
    module_logger.info('Computing the distance matrix using %i cores' % num_cores)
    dist_mat, shift_mat, _ = distance_matrix(
        stream_list = stream_list,
        shift_len = shift_len,
        cores=num_cores,
        replace_nan_distances_with=replace_nan_distances_with,
        allow_individual_trace_shifts=allow_individual_trace_shifts
    )
    # Handle nan entries in dist_mat
    dist_mat = handle_distmat_nans(dist_mat, replace_nan_distances_with=replace_nan_distances_with)
    return dist_mat, shift_mat

def cluster_correlated_templates(
        dist_mat, 
        template_list,
        corr_thresh=0.5,
        method='single', 
        metric='euclidian', 
        optimal_ordering=False,
        criterion='distance'):
    if not np.isfinite(dist_mat).all():
        module_logger.critical('dist_mat has non-finite entries - filling needs to be applied')
    if not 0 < corr_thresh < 1:
        module_logger.critical('corr_thresh is outside bounds of (0,1)')
    dist_vect = squareform(dist_mat)
    Z = linkage(dist_vect, method=method, metric=metric, optimal_ordering=optimal_ordering)
    indices = fcluster(Z, t=1 - corr_thresh, criterion=criterion)
    # Indices start at 1...
    group_ids = list(set(indices))  # Unique list of group ids
    module_logger.info(' '.join(['Found', str(len(group_ids)), 'groups']))
    # Convert to tuple of (group id, stream id)
    indices = [(indices[i], i) for i in range(len(indices))]
    # Sort by group id
    indices.sort(key=lambda tup: tup[0])
    groups = []
    module_logger.info('Extracting and grouping')
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
    return groups

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
        templates, show=False, corr_thresh=0.3, shift_len=1.,
        allow_individual_trace_shifts=False, fill_value=None,
        cores='all', save_path='.', save_subtribes=True,
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
    module_logger.info('Computing the distance matrix using %i cores' % num_cores)
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
        module_logger.info('Saved the distance matrix as dist_mat.npy')
        np.save(os.path.join(save_path, 'shift_mat.npy'), shift_mat)
        module_logger.info('Saved the shift matrix as shift_mat.npy')
    
    # Calculate linkage
    module_logger.info('Computing linkage')
    dist_mat = handle_distmat_nans(dist_mat,
                                   replace_nan_distances_with=fill_value)
    dist_vec = squareform(dist_mat)
    Z = linkage(dist_vec, **kwargs)

    # Get the indices of the groups
    module_logger.info('Clustering')
    indices = fcluster(Z, t=1 - corr_thresh, criterion='distance')
    # Indices start at 1...
    group_ids = list(set(indices))  # Unique list of group ids
    module_logger.info(' '.join(['Found', str(len(group_ids)), 'groups']))
    # Convert to tuple of (group id, stream id)
    indices = [(indices[i], i) for i in range(len(indices))]

    # Sort by group id
    indices.sort(key=lambda tup: tup[0])
    groups = []
    module_logger.info('Extracting and grouping')
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
        module_logger.info('Plotting the dendrogram')
        dendrogram(Z, color_threshold=1 - corr_thresh,
                   distance_sort='ascending')
        plt.show()
    return groups


def tribe_cluster(
        tribe, save_path='./cluster_results',
        corr_thresh=0.3, shift_len=1., fill_value = 1, 
        allow_individual_trace_shifts=False, cores='all', 
        show=True, method='linear',
        metric='euclidian', optimal_ordering=False):
    """Conduct waveform-based clustering of :class:`~eqcorrscan.core.match_filter.template.Template`
    objects using the :meth:`~eqcorrscan.util.clustering.cluster` method and record the parameter
    inputs and results (including the distance matrix) in a specified directory.

    We provide an extended set of input parameters to encompass patameters for this method,
    and the underlying :meth:`~eqcorrscan.util.clustering.cluster` and :meth:`~scipy.cluster.hierarchy.linkage`
    methods. 

    `tribe_cluster` Parameters
    --------------------------
    :param tribe: tribe containing templates
    :type: eqcorrscan.core.match_filter.tribe.Tribe
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
    :param method: method to use for linkage, defaults to 'linear'
    :type method: str, optional
    :param metric: metric to use for linkage, defaults to 'euclidian'
    :type metric: str, optional
    :param optimal_ordering: should optimal ordering be used? Defaults to False
    :type optimal_ordering: bool, optional

    Files Generated
    ---------------
     - **params.csv** -- all input parameters to this method as a CSV
        with parameter names in the left column and values in the right column
     - **groups_cct###.csv** -- names of templates (left column) and their group number (right column)
        based on the provided **corr_thresh** value (_cct### = _cct{corr_thresh: %.2f})
     - **dist_mat_{fill_value}.npy*** -- distance matrix saved by 

    Output
    ------
    :return:
     - **subtribes** (*dict*) - dictionary of eqcorrscan.Tribe objects that contain
        grouped templates, keyed to group number

    """    
    if not isinstance(tribe, Tribe):
        raise TypeError('input "tribe" must be type eqcorrscan.Tribe')
    elif len(tribe) < 2:
        raise ValueError(f'tribe of length {len(tribe)} cannot be clustered. 2 minimum.')
    else:
        pass

    # Compatability check for save_path
    if not isinstance(save_path, str):
        raise TypeError
    elif not os.path.exists(save_path):
        os.makedirs(save_path)
    else:
        pass

    # Compatability checks for fill_value
    if fill_value in ['mean','min']:
        pass
        dist_mat_name = f'dist_mat_{fill_value}.npy'
    elif isinstance(fill_value, float):
        pass
        dist_mat_name = 'dist_mat_%.2f.npy'%(fill_value)
    elif fill_value is None:
        dist_mat_name = 'dist_mat_None.npy'
    else:
        raise ValueError(f'fill_value {fill_value} not supported.')
    
    # Pull together `cluster` key-words
    kwargs = {'show': show,
              'corr_thresh': corr_thresh,
              'shift_len': shift_len,
              'allow_individual_trace_shift': allow_individual_trace_shifts,
              'save_corrmat': True,
              'replace_nan_distances_with': fill_value,
              'cores': cores,
              'method': method,
              'metric': metric,
              'optimal_ordering': optimal_ordering}
    
    # Save kwargs as parameters
    with open(os.path.join(save_path,'params.csv'), 'w') as _f:
        for _k, _v in kwargs.items():
            _f.write(f'{_k},{_v}\n')
    
    # Group Templates
    template_list = _compose_template_list(tribe)
    
    # RUN CLUSTERING
    groups = cluster(template_list, **kwargs)
    
    # Put templates into groups
    subtribes = {}
    for _e, group in enumerate(groups):
        itribe = Tribe()
        for entry in group:
            name = entry[1]
            itemplate = tribe.select(name)
            itribe += itemplate
        subtribes.update({_e: itribe})

    # Move dist_mat.npy & re-label
    src = os.path.join(Path().cwd(), 'dist_mat.npy')
    dest = os.path.join(save_path,dist_mat_name)
    if os.path.isfile(src):
        os.rename(src, dest)

    # Save group membership as CSV
    with open(os.path.join(save_path, f'groups_cct{corr_thresh:.2f}'), 'w') as _f:
        _f.write('template_name,group_id\n')
        for _e, group in enumerate(groups):
            for entry in group:
                _f.write(f'{entry[1]},{_e}\n')
            
    return subtribes


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



