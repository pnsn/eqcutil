
import os, logging, tarfile, shutil, pickle, tempfile, glob
from pathlib import Path

import numpy as np
import pandas as pd

from obspy import read_events, read
from obspy.core.event import Catalog
from eqcorrscan import Tribe, Template
import eqcorrscan.utils.clustering as euc
from eqcorrscan.core.match_filter.helpers import _safemembers, _par_read

from eqcorrscan_utils.augment.template import rename_templates, deduplicate_names
from eqcorrscan_utils.io.template_cluster import save_cluster_results
from eqcorrscan_utils.util.decorators import catch_kwargs
from eqcorrscan_utils.visualize.snuffler import plant

Logger = logging.getLogger(__name__)


# def corr_cluster(self, corr_thresh=0.3, shift_len=1., allow_individual_trace)

class ClusteringTribe(Tribe):
    """An augmentation of the :class:`~eqcorrscan.Tribe` class that
    extends EQcorrscan clustering and stacking class-methods beyond
    those in the current distribution of EQcorrscan

    """    
    def __init__(self, templates=[], use_evid_names=False, **options):
        if not hasattr(Tribe,'snuffle'):
            plant()
        if isinstance(templates, Template):
            templates = [templates]
        elif isinstance(templates, list):
            if all(isinstance(_t, Template) for _t in templates):
                pass
            else:
                raise TypeError
        if use_evid_names:
            if 'include_contributor_code' in options.keys():
                templates = rename_templates(include_contributor_code=options['include_contributor_code'])
                Logger.info('renamed templates using event resource_id elements')
            else:
                templates = rename_templates(templates)
    
        if len(templates) > {_t.name for _t in templates}:
            templates = deduplicate_names(templates)

        super().__init__(templates=templates)

        self.clusters = {}
        self.stacks = {}
        self.cluster_kwargs = {}
        self.stack_kwargs = {}

    def get_summary(self):
        index = [_t.name for _t in self]
        columns = list(self.clusters.key())
        if len(columns) > 0:
            data = np.full(shape=(len(index), len(columns)), fill_value=np.nan)
        df = pd.DataFrame(data=data, index=index, columns=columns)
        # Iterate over cluster type
        for ctype, cgroups in self.clusters.items():
            # Iterate over subtribe number and subtribe
            for _stn, tribe in cgroups.items():
                # Iterate over template in tribe
                for template in tribe:
                    # Write subgroup number to name,cluster_type position
                    df.loc[template.name, ctype] = _stn
        return df

    summary = property(get_summary)

    def get_template_list(self, use_name=False):
        """Create a template_list input for 
        :meth:`~eqcorrscan.utils.clusering.cluster`
        from this ClusteringTribe with the option
        to use template names instead of position
        values.

        :return: _description_
        :rtype: _type_
        """
        if use_name:
            return [(_t.st, _t.name) for _t in self]      
        else:       
            return [(_t.st, _e) for _e, _t in enumerate(self)]
    
    template_list = property(get_template_list)


    def cluster(self, method, **kwargs):
        """Extended wrapper for EQcorrscan Template correlation methods
        In addition to the original options of clustering using catalog
        methods (space_cluster and (space_time_cluster) this method now
        also permits use of :meth:`~

        :param method: _description_
        :type method: _type_
        """        
        if method in ['space_cluster','space_time_cluster']:
            tribes = Tribe.cluster(self, method, **kwargs)
        elif method == 'correlation_cluster':
            groups = euc.cluster(self.get_template_list(), **kwargs)
            tribes = []
            for group in groups:
                templates = []
                for entry in group:
                    templates.append(self.templates[entry[1]])
                tribes.append(templates)

        self.cluster_kwargs.update({method: kwargs})
        self.clusters.update({method: {_e: self.__class__(templates=_t) for _e, _t in enumerate(tribes)}})
    
    def write(self, filename, compress=True, catalog_format="QUAKEML"):
        """
        Write the clusteringtribe to a file using tar archive formatting.

        :type filename: str
        :param filename:
            Filename to write to, if it exists it will be appended to.
        :type compress: bool
        :param compress:
            Whether to compress the tar archive or not, if False then will
            just be files in a folder.
        :type catalog_format: str
        :param catalog_format:
            What format to write the detection-catalog with. Only Nordic,
            SC3ML, QUAKEML are supported. Note that not all information is
            written for all formats (QUAKEML is the most complete, but is
            slow for IO).

        .. rubric:: Example

        >>> tribe = Tribe(templates=[Template(name='c', st=read())])
        >>> tribe.write('test_tribe')
        Tribe of 1 templates
        >>> tribe.write(
        ...    "this_wont_work.bob",
        ...    catalog_format="BOB") # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
        TypeError: BOB is not supported
        """
        from eqcorrscan.core.match_filter import CAT_EXT_MAP

        if catalog_format not in CAT_EXT_MAP.keys():
            raise TypeError("{0} is not supported".format(catalog_format))
        dirname, ext = os.path.splitext(filename)

        # Make directory if it doesn't exist
        if not os.path.isdir(dirname):
            os.makedirs(dirname)
        self._par_write(dirname)

        # Compose tribe catalog
        tribe_cat = Catalog()
        for t in self.templates:
            if t.event is not None:
                # Check that the name in the comment matches the template name
                for comment in t.event.comments:
                    if not comment.text:
                        comment.text = "eqcorrscan_template_{0}".format(t.name)
                    elif comment.text.startswith("eqcorrscan_template_"):
                        comment.text = "eqcorrscan_template_{0}".format(t.name)
                tribe_cat.append(t.event)
        # Write catalog to disk
        if len(tribe_cat) > 0:
            tribe_cat.write(
                os.path.join(dirname, 'tribe_cat.{0}'.format(
                    CAT_EXT_MAP[catalog_format])), format=catalog_format)
        # Write template streams to disk
        for template in self.templates:
            template.st.write(
                os.path.join(dirname, '{0}.ms'.format(template.name)),
                format='MSEED')
        # ADDED BY NTS - write summary to disk
        self.summary.to_csv(os.path.join(dirname,'cluster_summary.csv'), header=True, index=True)
        # Write clustering params to CSV
        for _k, _v in self.cluster_kwargs:
            with open(os.path.join(dirname, f'{_k}_kwargs.csv'), 'w') as file:
                for _l, _w in _v:
                    if isinstance(_w, str):
                        file.write(f'{_l},{_w}\n')
                    else:
                        file.write(f'{_l},{repr(_w)}\n')
        if compress:
            if not filename.endswith(".tgz"):
                Logger.info("Appending '.tgz' to filename.")
                filename += ".tgz"
            with tarfile.open(filename, "w:gz") as tar:
                tar.add(dirname, arcname=os.path.basename(dirname))
            shutil.rmtree(dirname)
        return self

    def read(self, filename):
        """
        Read a clustertribe of templates from a tar formatted file.

        :type filename: str
        :param filename: File to read templates from.

        .. rubric:: Example

        >>> tribe = Tribe(templates=[Template(name='c', st=read())])
        >>> tribe.write('test_tribe')
        Tribe of 1 templates
        >>> tribe_back = Tribe().read('test_tribe.tgz')
        >>> tribe_back == tribe
        True
        >>> # This can also read pickled templates
        >>> import pickle
        >>> with open("test_tribe.pkl", "wb") as f:
        ...    pickle.dump(tribe, f)
        >>> tribe_back = Tribe().read("test_tribe.pkl")
        >>> tribe_back == tribe
        True
        """
        if filename.endswith(".pkl"):
            with open(filename, "rb") as f:
                self.__iadd__(pickle.load(f))
            return self
        with tarfile.open(filename, "r:*") as arc:
            temp_dir = tempfile.mkdtemp()
            arc.extractall(path=temp_dir, members=_safemembers(arc))
            tribe_dir = glob.glob(temp_dir + os.sep + '*')[0]
            self._read_from_folder(dirname=tribe_dir)
        shutil.rmtree(temp_dir)
        # Assign unique ids
        self.__unique_ids()
        return self

    def _read_from_folder(self, dirname):
        """
        Internal folder reader.

        :type dirname: str
        :param dirname: Folder to read from.
        """
        templates = _par_read(dirname=dirname, compressed=False)
        t_files = glob.glob(dirname + os.sep + '*.ms')
        tribe_cat_file = glob.glob(os.path.join(dirname, "tribe_cat.*"))
        cluster_file = glob.glob(os.path.join(dirname,'cluster_summary.csv'))
        cluster_kwarg_files = glob.glob(os.path.join(dirname,'*_kwargs.csv'))
        if len(tribe_cat_file) != 0:
            tribe_cat = read_events(tribe_cat_file[0])
        else:
            tribe_cat = Catalog()
        previous_template_names = [t.name for t in self.templates]
        for template in templates:
            if template.name in previous_template_names:
                # Don't read in for templates that we already have.
                continue
            for event in tribe_cat:
                for comment in event.comments:
                    if comment.text == 'eqcorrscan_template_' + template.name:
                        template.event = event
            t_file = [t for t in t_files
                      if t.split(os.sep)[-1] == template.name + '.ms']
            if len(t_file) == 0:
                Logger.error('No waveform for template: ' + template.name)
                continue
            elif len(t_file) > 1:
                Logger.warning('Multiple waveforms found, using: ' + t_file[0])
            template.st = read(t_file[0])
        # Remove templates that do not have streams
        templates = [t for t in templates if t.st is not None]
        self.templates.extend(templates)
        # Re-constitute groups
        summary = pd.read_csv(cluster_file, index_col=[0])
        for name, row in summary.iterrows():
            for col in row.columns:
                if col not in self.clusters.keys():
                    self.clusters.update({col: self.__class__()})
                    self.clusters[col] += self.select(name)
        # Reconstitute processing information
        for ckf in cluster_kwarg_files:
            path, name = os.path.split(ckf)
            name, ext = os.path.splitext(name)
            ctype = name[:-7]
            self.cluster_kwargs.update({ctype: {}})
            df = pd.read_csv(ckf, index_col=[0])
            for _k, _r in df.iterrows():
                self.cluster_kwargs[ctype].update({_k, _r.values[0]})
                
        return

    # def corr_cluster(self, savedir='cluster_results', show=False, corr_thresh=0.3,
    #             shift_len=1., allow_individual_trace_shifts=False, dist_nan_fill=None,
    #            cores='all', method='linear', metric='euclidian',
    #            optimal_ordering=False, save_subtribes=False):
    #     """Wraps the :meth:`~eqcorrscan.util.clustering.cluster` method and 
    #     added saving and formatting functionalities provided by the eqcorrscan_utils
    #     package. 

    #     :param savedir: where clustering results are saved, defaults to 'cluster_results'
    #     :type savedir: str, optional
    #     :param show: should the sub-call of `cluster` display the groups as a dendrogram? Defaults to False
    #     :type show: bool, optional
    #     :param corr_thresh: correlation threshold for grouping, defaults to 0.3
    #     :type corr_thresh: float, optional
    #     :param shift_len: number of seconds templates are allowed to be shifted
    #         durring cross correlations, defaults to 1.
    #     :type shift_len: float, optional
    #     :param allow_individual_trace_shifts: allow each trace in a template to be
    #         shifted during cross correlations? Defaults to False
    #     :type allow_individual_trace_shifts: bool, optional
    #     :param dist_nan_fill: fill value for NaN entries in the distance matrix
    #         calculated by `cluster` (its replace_nan_distances_with arg),
    #         defaults to None.
    #         Supported Values:
    #          - 'mean'
    #          - 'min'
    #          - [0, 1]
    #          also see :meth:`~eqcorrscan.utils.clustering.cluster`
    #     :type dist_nan_fill: str, float, or NoneType, optional
    #     :param cores: number of cores to use for `cluster`, defaults to 'all'
    #     :type cores: str, optional
    #     :param method: linkage method, defaults to 'linear'
    #         also see :meth:`~scipy.clustering.hierarchy.linkage
    #     :type method: str, optional
    #     :param metric: linkage metric, defaults to 'euclidian'
    #     :type metric: str, optional
    #     :param optimal_ordering: linkage optimal_ordering, defaults to False
    #     :type optimal_ordering: bool, optional
    #     :param save_subtribes: should subtribes be saved in savedir?, defaults to False
    #     :type save_subtribes: bool, optional
    #     """



    #     # Create Template List
    #     template_list = [(_t.st, _e) for _e, _t in enumerate(self)]
    #     # Run Clustering
    #     groups = euc.cluster(
    #         template_list,
    #         show=show,
    #         corr_thresh=corr_thresh,
    #         shift_len=shift_len,
    #         allow_individual_trace_shifts=allow_individual_trace_shifts,
    #         replace_nan_distances_with=dist_nan_fill,
    #         cores=cores,
    #         method=method,
    #         metric=metric,
    #         optimal_ordering=optimal_ordering)
    #     # Populate/re-initialize subtribes
    #     self.clusters = {}
    #     if savedir:
    #         subtribes = save_cluster_results(self, groups, savedir=savedir, save_subtribes=save_subtribes)
    #         for _k, _v in subtribes.items():
    #             self.clusters.update({_k, self.__class__(templates=_v)})
    #     else:
    #         for _e, group in groups:
    #             self.clusters.update({_e: self.__class__()})
    #             for entry in group:
    #                 self.clusters[_e] += self.templates[entry[1]]
    #     # Load the distance matrix into this 
    #     self.dist_mat = np.load(Path().cwd() / 'dist_mat.npy')
    #     self.clustering_threshold=corr_thresh
    #     self.dist_nan_fill = dist_nan_fill
    #     self.savedir = savedir
    

    # def write(self, savedir):
    #     if self.groups is None:
    #         Tribe.write(self, savedir, 'templates.tgz')
    #     gfile = os.path.join(savedir, 'groups.csv')
        


    # def load(self, loaddir):
    #     gfile = os.path.join(loaddir, 'groups.csv')
    #     pfile = os.path.join()
    #     if not os.path.isfile(os.path.join(loaddir,'groups.csv')):
    #         Logger.critical(f'could not find groups.csv in {loaddir}')
    #     else:
    #         with open()
