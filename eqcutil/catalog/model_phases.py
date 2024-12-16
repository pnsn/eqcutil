"""
:module: eqcorrscan_utils.cake
:auth: Nathan T. Stevens
:email: ntsteven@uw.edu
:org: Pacific Northwest Seismic Network
:license: GNU GPLv3
:purpose:
    This module provides a more-approachable implementation of the pyrocko.cake
    module for ray-tracing in a 1-D velocity structure and integration with
    ObsPy seismic event and station metadata objects that are used with EQcorrscan

TODO: Consider making this a class that contains
"""
import logging
import pandas as pd
from obspy import Inventory
from obspy.core.event import Pick, Arrival, WaveformStreamID
from obspy.geodetics import locations2degrees
from pyrocko import cake

def model_arrivals(origin, inventory, model_name='P4', phases=['P','S']):
    """Model raypaths and travel-times between a seismic event
    and a set of seismic receivers for a prescribed PNSN mdoel
    and phase types.
    
    Parameters
    ----------
    :param origin: seismic event object
    :type origin: obspy.core.event.Origin
    :param inventory: seismic station inventory containing station objects
    :type inventory: obspy.core.inventory.Inventory
    :param model_name: name of the PNSN model to use, defaults to 'P4'
    :type model_name: str, optional
    :param phases: name of phases to use, defaults to ['P','S']
    :type phases: list, optional
    :return:
     - **summary** (*pandas.DataFrame*) -- dataframe summarizing the
        modeled ray-paths and parameters needed for static corrections
    
    Columns
    -------
     - station: station code
     - sta dz m: Station elevation offset in meters (0 - elevation)
     - orig dz m: Origin depth offset in meters (positive is above sea-level)
     - ray parameter s/deg: ray parameter in seconds per degree
     - offset km: epicenter-receiver offset in kilometers
     - takeoff angle: takeoff angle in degrees
     - incidence angle: incidence angle in degrees
     - transmission efficiency: percent efficiency of transmission due to path effects
     - spreading efficiency: percent efficincy of transmission due to spreading effects
     - path description: description of the ray-path:
            Given Phase    (Used Phase)   (start_layer - turning_layer - ending_layer)
  
    """    
    model = make_model(name=model_name)
    raypaths, e_offset = model_raypaths(model, origin, inventory, phases=phases)
    summary = ray_summary(raypaths, origin, inventory, e_offset=e_offset)
    return summary

## MODEL CONSTRUCTION METHODS ##

def get_pnsn_model(name='P4'):
    """Load 1-D velocity Models used by the PNSN

    :param name: velocity model name, defaults to 'P4'
    :type name: str, optional
        Supported Values:
         - 'P4'/'P5': Puget Sound / Western Washington Region (updated Sept 2020)
         - 'C4': Cascades Region (updated Sept 2020)
    :raises NotImplementedError: _description_
    :return: _description_
    :rtype: _type_
    """    
    if name in ['P4', 'P5']:
        Vp = [5.4, 6.38, 6.59, 6.73, 6.86, 6.95, 7.8]
        VpVs = 1.78
        Ztop = [0., 4., 9., 16., 20., 25., 41.]
    elif name == 'C4':
        Vp = [5.1, 6.0, 6.6, 6.8, 7.1, 7.8]
        VpVs = 1.78
        Ztop = [0., 1., 10., 18., 34., 43.]
    else:
        raise NotImplementedError(f'name "{name}" not supported. Options: P5, C4')
    return Vp, VpVs, Ztop

def create_1d_model(Vp, VpVs, Ztop, padding_meters=1e5):
    """Create a 1D layered velocity model with homogeneous velocity
    layers (i.e., velocity gradient is 0)
    
    :class:`~pyrocko.cake.LayeredModel` object containing
    :class:`~pyrocko.cake.HomogeneouseLayer` elements.

    :param Vp: layer P-wave velocities in km/sec
    :type Vp: list of float-like
    :param VpVs: VpVs ratio contstant for this model
    :type VpVs: float-like
    :param Ztop: layer top depths in km
    :type Ztop: list of float-like
    :param padding_meters: thickness of the final layer in meters
        defaults to 100000 (100 km).
    :type padding_meters: float-like
    :return: 
     - **model** (*pyrocko.cake.LayeredModel*) - model object
    """    
    model = cake.LayeredModel()
    if len(Vp) != len(Ztop):
        raise ValueError(f'{len(Vp)} != {len(Ztop)}')
    else:
        pass

    for _e in range(len(Vp)):
        vp = Vp[_e]*1e3
        vs = vp/VpVs
        material = cake.Material(vp=vp, vs=vs)
        ztop = Ztop[_e]*1e3
        if _e + 1 < len(Vp):
            zbot = Ztop[_e + 1]*1e3
        else:
            zbot = Ztop[_e]*1e3 + padding_meters
        layer = cake.HomogeneousLayer(ztop=ztop, zbot=zbot, m=material)
        model.append(layer)

    return model

def make_model(name='P4'):
    """Wrapper script for constructing a LayeredModel
    for a designated 

    :param name: name of the PNSN model to construct, defaults to 'P4'
    :type name: str, optional
    :return: 
        - **model** (*pyrocko.cake.LayeredModel*) - model object
    """    
    Vp, VpVs, Ztop = get_pnsn_model(name=name)
    return create_1d_model(Vp, VpVs, Ztop)



#############################
## RAYPATH MODELING METHOD ##
#############################

def model_raypaths(model, origin, inventory, phases=['P','S']):
    results = {}
    # Get Origin Hypocentral parameters
    olat = origin.latitude
    olon = origin.longitude
    odep = origin.depth
    oele = -1.*origin.depth
    time = origin.time
    # If origin is below sea-level, set first guess at e_offset as sea-level
    if oele < 0:
        e_offset = 0.
    # If origin is above sea-level, take that as the first guess at e_offset
    else:
        e_offset = oele # [m] above sea-level
    distances = set([])
    # Iterate across networks & stations
    for net in inventory.networks:
        for sta in net.stations:
            slat = sta.latitude
            slon = sta.longitude
            sele = sta.elevation
            # If station sits above top model, adjust e_offset
            if sele > e_offset:
                e_offset = sele
            # Calculate lateral offset
            delta = locations2degrees(olat, olon, slat, slon)
            distances.add(delta)
    # Model offsets
    arrivals = model.arrivals(distances=list(distances),
                              zstart=odep - e_offset,
                              zstop=0.)
    # Zip modeling results
    ddict = dict(zip(distances, arrivals))
    # Associate stations & results
    for net in inventory.networks:
        for sta in net.stations:
            slat = sta.latitude
            slon = sta.longitude
            sele = sta.elevation
            # If station sits above top model, adjust e_offset
            if sele > e_offset:
                e_offset = sele
            # Calculate lateral offset
            delta = locations2degrees(olat, olon, slat, slon)
            results.update({sta.code: ddict[delta]})

    return (results, e_offset)



def model_raypaths_simple(model, origin, inventory, phases=['P','S']):
    """
    Model the P-wave travel times at stations in the provided inventory
    from an earthquake origin. This model assumes that stations are at 
    sea-level and any sources above sea-level (i.e. surface
    events) are pinned to sea-level (depth = 0) 

    Additional station corrections should be applied to modeled travel times
    to correct for 

    :param model: 1-D velocity model
    :type model: pyrocko.cake.LayeredModel
    :param origin: origin solution object
    :type origin: obspy.core.event.Origin
    :param inventory: station inventory
    :type inventory: obspy.core.inventory.Inventory

    :returns:
     - **results** (*dict*) -- dictionary keyed by station code (Net.Sta)
        with list values comprising sets of :class:`~pyrocko.cake.Ray` objects
    """    
    results = {}
    # Iterate across networks in inventory
    for network in inventory:
        # Iterate across stations
        for station in network.stations:
            # Get station-epicenter distance
            dist_deg = locations2degrees(
                origin.latitude,
                origin.longitude,
                station.latitude,
                station.longitude
            )
            if origin.depth > 0:
                zstart = origin.depth
            else:
                zstart = 0

            arrivals = model.arrivals(distances=[dist_deg],
                                    zstart=zstart,
                                    zstop=0.,
                                    phases=phases)
            if arrivals:
                results[station.code] = arrivals
    return results

########################
## CONVERSION METHODS ##
########################
def make_wfid(nslc, resource_uri=None):
    """Convert a SEED channel code string into a 
    :class:`~obspy.core.event.WaveformStreamID`
    object.

    :param nslc: SEED channel code
    :type nslc: str
    :return: 
     - **waveform_id** (*obspy.core.event.WaveformStreamID*) -- waveform stream ID object
    """    
    return WaveformStreamID(seed_string=nslc, resource_uri=resource_uri)

def ray2pick(ray, wfid, origin):
    """
    Create a :class:`~obspy.core.event.Pick` representation of
      a :class:`~pyrocko.cake.RayPath` object
    
    :param ray: raypath object
    :type ray: pyrocko.cake.RayPath
    :param wfid: waveform stream ID
    :type wfid: obspy.core.event.WaveformStreamID
    :param origin: event origin used to generate **ray**
    :type origin: obspy.core.event.Origin

    """    
    t0 = origin.time
    ta = t0 + ray.t
    pick = Pick(time=ta,
                phase_hint=ray.phase.path.given_name(),
                waveform_id=wfid,
                evaluation_mode='automatic',
                time_errors=None)
    return pick


def ray_summary(results, origin, inventory):
    """Create a human-readable synopsis of outputs from 
    :meth:`~.model_raypaths` in the form of a :class:`~pandas.DataFrame`
    object.

    Columns:
     - station: station code
     - sta dz m: Station elevation offset in meters (0 - elevation)
     - orig dz m: Origin depth offset in meters (positive is above sea-level)
     - ray parameter s/deg: ray parameter in seconds per degree
     - offset km: epicenter-receiver offset in kilometers
     - takeoff angle: takeoff angle in degrees
     - incidence angle: incidence angle in degrees
     - transmission efficiency: percent efficiency of transmission due to path effects
     - spreading efficiency: percent efficincy of transmission due to spreading effects
     - path description: description of the ray-path:
            Given Phase    (Used Phase)   (start_layer - turning_layer - ending_layer)

    :param results: modeled arrivals for one or more stations
    :type results: dict
    :param origin: Origin object used to model ray-paths
    :type origin: obspy.core.event.Origin
    :param inventory: station inventory used to model ray-paths
    :type inventory: obspy.core.inventory.Inventory

    :return: 
        - **df** (*pandas.DataFrame*) summary dataframe
    """    
    holder = []
    if origin.depth < 0:
        dz_orig = 0 - origin.depth
    else:
        dz_orig = origin.depth

    for _sta, _arrivals in results.items():
        iinv = inventory.select(station=_sta)
        sta = iinv[0].stations[0] 
        dz_sta = 0 - sta.elevation
        for _arr in _arrivals:
            line = [_sta,
                    dz_sta,
                    dz_orig,
                    _arr.path.phase.given_name(),
                    _arr.p/cake.r2d,
                    _arr.x*(cake.d2r*cake.earthradius/cake.km),
                    _arr.t,
                    _arr.takeoff_angle(),
                    _arr.incidence_angle(),
                    100.*_arr.efficiency(),
                    100.*_arr.spreading()*_arr.surface_sphere(),
                    _arr.path.__str__(p=_arr.p)]
            holder.append(line)
            

    df = pd.DataFrame(data=holder,
                      columns=['station','sta dz m','orig dz m','phase','ray parameter s/deg','offset km','travel time sec','takeoff angle','incidence angle',
                               'transmission efficiency','spreading efficiency','path description'])
    
    df.sort_values(by=['offset km','travel time sec'], ascending=True, inplace=True)

    return df