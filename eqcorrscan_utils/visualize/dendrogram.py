
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform

from eqcorrscan.utils.clustering import cluster, handle_distmat_nans, distance_matrix, fcluster



def reconstitute_linkage(loaddir, **kwargs):
    """Re-load clustering data saved by :meth:`~.save_cluster_results`
    and reconstitue the linkage that serves as an input for creating
    dendrograms 
    
    
    a saved `dist_mat.npy`-type file saved by
    :meth:`~eqcorrscan.util.cluster.clustering` and re-create
    the `Z` input to :meth:`~scipy.cluster.dendrogram` via the
    :meth:`~scipy.cluster.hierarchical.linkage` method and provide
    an option for a NaN fill_value

    **fill_value** corresponds to the `replace_nan_distances_with` parameter
    in :meth:`~eqcorrscan.util.cluster.clustering`.
    
    .. rubric:: Notes on `fill_value`
        Explainer for fill_value selection:
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
    if os.path.isfile(os.path.join(loaddir, 'dist_mat.npy')):
        dist_mat = np.load(os.path.join(loaddir, 'dist_mat.npy'))
    else:
        Logger.critical(f'dist_mat.npy not found in {loaddir}')
    if os.path.isfile(os.path.join(loaddir, 'cluster_params.csv')):
        

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