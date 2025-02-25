from collections import defaultdict
import numpy as np
from obspy import Catalog, Trace, Stream, UTCDateTime
from obspy.geodetics import locations2degrees, degrees2kilometers
from eqcorrscan import Tribe, Template


class Grouper(object):
    """Rewrite of the EQcorrscan :class:`~eqcorrscan.Tribe` class intended
    to better facilitate class-method based clustering/grouping routines
    

    :param object: _description_
    :type object: _type_
    """    
    def __init__(self, wavebank_basepath, eventbank_basepath):
        self.catalog = defaultdict(Catalog)
        self.channels = defaultdict(Stream)
        self.index ={'name': defaultdict(dict), 
                     'trace': defaultdict(dict)}
        # Arrays (R3) - evid x evid x trace_id
        self.dist = np.array([[[]]])
        self.dt = np.array([[[]]])
        self.corr = np.array([[[]]])
        self.shift = np.array([[[]]])
        self._has_proc_info = False
        for _t in templates:
            self._add_template(_t)

    def _add_template(self, other, **options):
        # Ensure other is eqcorrscan Template
        if not isinstance(other, Template):
            raise TypeError
        # Get name
        name = other.name
        # If new event
        if name not in self.index['name'].keys():
            self.index['name'][name] = len(self.index['name'])
        # Get ij coordinate
        ij = self.index['name'][name]
        # Iterate across traces in other
        for tr in other.st:
            # If new trace ID, add new index entry
            if tr.id not in self.index['trace'].keys():
                self.index['trace'][tr.id] = len(self.index['trace'])
            kk = self.index['trace'][tr.id]
            new_dist = np.full(shape=(ij,ij,kk), fill_value=np.nan)
            new_dt = np.full(shape=(ij,ij,kk), fill_value=np.nan)
            new_

            # Iterate across events
            for _ij in self.index['name'].values():
                np.pad(self.dist, ((0,1),(0,1)),
                       mode='constant',
                       constant_values=np.nan)
                np.pad(self.dt, ((0,1),(0,1)),
                       mode='constant',
                       constant_values=np.nan)
                np.pad(self.corr, ((0,1),(0,1)),
                       mode='constant',
                       constant_values=np.nan)
                np.pad(self.shift, ((0,1),(0,1)),
                       mode='constant',
                       constant_values=np.nan)

                

        

        # Ensure processing info matches / ingest processing info
        if not self._has_proc_info:
            self.proc_info = {_k: getattr(other, _k) for _k in ['lowcut','highcut','samp_rate','filt_order','process_length','prepick']}
            self._has_proc_info = True
        else:
            for _k, _v in self.proc_info.items():
                if _v != getattr(other, _k):
                    raise AttributeError(f'{_k} has mismatch')
        # Capture/Sort Traces, EVIDs, and Trace IDs
        self.index[other.name]
        if other.name in self.index.keys():
            for tr in other.st:
                if tr.id in self.index[other.name].keys():

        # Is event ID in 
        if len(self.index) == 0:
            for 
            self.index[other.name]

        
        if other.name not in self.events.keys():
            self.events[other.name] += other.event
            self.streams[other.name] += other.st
            new_event = True
        else:
            new_event = False

        if new_event:
            pass
        else:
            for tr in 

        default_kwargs = {
            'shift_len': other.prepick*2,
            'allow_individual_trace_shifts': True,
            'xcorr_func': 'fftw',
            'cores': 1
            }
        default_kwargs.update(options)
        scale = len(self)
        for key in ['dist','dt','corr']:
            if scale == 0:
                self.__setattr__(key, np.zeros(shape=(1,1)))
            elif:
                _array = self.__getattr__(key)
                # Expand array
                _array = np.pad(_array, ((0, 1), (0, 1)), mode='constant', constant_values=np.nan)
                for ii in range(scale):
                    itemp = self.templates[ii]
                    iprefor = itemp.event.preferred_origin()
                    oprefor = other.event.preferred_origin()
                    if key == 'dist':
                        # Get degrees distance
                        _d = locations2degrees(iprefor.latitude,
                                               iprefor.longitude,
                                               oprefor.latitude,
                                               oprefor.longitude)
                        # Convert to KM
                        _d = degrees2kilometers(_d)
                        # Incorporate depth difference
                        dz = iprefor.depth - oprefor.depth
                        _d = np.sqrt(_d**2 + dz**2)
                    elif key == 'dt':
                        _d = np.abs(iprefor.time - oprefor.time)
                    elif key == 'corr':
                        dist_list, shift_list = cross_chan_correlation(
                            other.st,
                            self.get_stream_list(),
                            **default_kwargs)
                    _array[ii, -1] = _d
                    _array[-1, ii] = _d
                
                _array[ii + 1, ii + 1] = 0.

                        
                


    def _calc_new_pairing(self, other):
        for key in ['dist','dt','corr']:

            if len(self) == 0:
                self.__setattr__(key, pd.DataFrame(data=[0], index=[other.name], columns=[other.name]))
            else:
                scale = len(self)
                _df = getattr(self, key)

        if len(self) == 0:
            for key in ['dist','dt','corr']:

            self.dist = pd.DataFrame()


def get_origin_comp(self, param='loc'):
        if param.lower() in ['distance','location','origin location','origin_location']:
            param='distance'
        elif param.lower() in ['time','origin time','origin_time']:
            param='time'
        elif param.lower() in ['origin']:
            raise NotImplementedError('TODO: implement space_time_cluster')
        else:
            raise ValueError(f'{param} not supported')
        # TODO: add functionality 
        locs = self.clusters.id_no.values
        # Get necessary scale for array
        scale = np.max(locs)
        # Add 1 for transition to length
        new_array = np.full(shape=(scale + 1, scale + 1), fill_value=np.nan)
        for _e, etemp in enumerate(self):
            for _f, ftemp in enumerate(self):
                ii = locs[_e]
                jj = locs[_f]
                if _e == _f:
                    new_array[ii, jj] = 0.
                elif _e < _f:
                    # Calculate distance in degrees
                    oe = etemp.event.preferred_origin()
                    of = ftemp.event.preferred_origin()
                    if param == 'distance':
                        ddeg = locations2degrees(oe.latitude, oe.longitude, of.latitude, of.longitude)
                        # Convert to kilometers
                        _d = degrees2kilometers(ddeg)
                    elif param == 'time':
                        _d = np.abs(oe.time - of.time)
                    else:
                        raise NotImplementedError
                    # Populate off-diagonal 
                    new_array[ii, jj] = _d
                    new_array[jj, ii] = _d
        if param == 'distance':
            self.dist_mat = new_array
        elif param == 'time':
            self.time_mat = new_array

    def calc_difftimes(self):
        self.get_origin_comp(param='time')
    
    def calc_distances(self):
        self.get_origin_comp(param='distance')
