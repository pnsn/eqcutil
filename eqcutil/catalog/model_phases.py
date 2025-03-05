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
from obspy.core.event import Pick, Arrival, WaveformStreamID, Origin, ResourceIdentifier
from obspy.geodetics import locations2degrees
from pyrocko import cake

# def model_arrivals(origin, inventory, model_name='P4', phases=['P','S']):
#     """Model raypaths and travel-times between a seismic event
#     and a set of seismic receivers for a prescribed PNSN mdoel
#     and phase types.
    
#     Parameters
#     ----------
#     :param origin: seismic event object
#     :type origin: obspy.core.event.Origin
#     :param inventory: seismic station inventory containing station objects
#     :type inventory: obspy.core.inventory.Inventory
#     :param model_name: name of the PNSN model to use, defaults to 'P4'
#     :type model_name: str, optional
#     :param phases: name of phases to use, defaults to ['P','S']
#     :type phases: list, optional
#     :return:
#      - **summary** (*pandas.DataFrame*) -- dataframe summarizing the
#         modeled ray-paths and parameters needed for static corrections
    
#     Columns
#     -------
#      - station: station code
#      - sta dz m: Station elevation offset in meters (0 - elevation)
#      - orig dz m: Origin depth offset in meters (positive is above sea-level)
#      - ray parameter s/deg: ray parameter in seconds per degree
#      - offset km: epicenter-receiver offset in kilometers
#      - takeoff angle: takeoff angle in degrees
#      - incidence angle: incidence angle in degrees
#      - transmission efficiency: percent efficiency of transmission due to path effects
#      - spreading efficiency: percent efficincy of transmission due to spreading effects
#      - path description: description of the ray-path:
#             Given Phase    (Used Phase)   (start_layer - turning_layer - ending_layer)
  
#     """    
#     model = make_model(name=model_name)
#     raypaths, e_offset = model_raypaths(model, origin, inventory, phases=phases)
#     summary = ray_summary(raypaths, origin, inventory)
#     return summary

def model_picks(origin, inventory, model_name='P4', phases=['P','S']):
    if not isinstance(origin, Origin):
        raise TypeError('origin must be type obspy.core.event.Origin')
    if not isinstance(inventory, Inventory):
        raise TypeError
    elif len(inventory) == 0:
        raise ValueError('Empty Inventory')
    if model_name not in ['P4']:
        raise ValueError(f'model_name "{model_name}" not supported')
    if isinstance(phases, str):
        phases = [phases]
    elif isinstance(phases, list):
        if all(isinstance(_e, str) for _e in phases):
            pass
        else:
            raise ValueError
    else:
        raise TypeError
    # Initialize velocity model object
    model = make_model(name=model_name)
    # Create holder for picks
    picks = []
    # Iterate across NSLC codes
    for nslc in inventory.get_contents()['channels']:
        # Get subset inventory for this specific stachan
        sta = nslc.split('.')[1]
        _inv = inventory.select(station=sta)
        # Model raypaths for this station
        results = model_raypaths(model, origin, _inv, phases=phases)
        # If modelling does not produce arrivals
        if results == []:
            breakpoint()
        # Convert rays to picks
        for net in _inv.networks:
            for sta in net.stations:
                wfid = WaveformStreamID(
                    network_code=net.code,
                    station_code=sta.code,
                    location_code=sta.channels[0].location_code,
                    channel_code=sta.channels[0].code
                )
                method_str = f'smi:local/pyrocko/cake/{model_name}'
                rays, e_offset, routing = results[f'{net.code}.{sta.code}']
                method_str = f'{method_str}/{e_offset:.3f}/{routing}'
                for ray in rays:
                    pick = ray2pick(ray, wfid, origin, method_str)
                    picks.append(pick)
        
    return picks

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
    """Updated version of model_raypaths that uses the Hypo2000 elevation adjustment
    approach to sources and/or events above sea-level. All receiver(s) in **inventory**
    and the input **origin** elevations are assessed for positive elevations / negative
    depths (respectively) and all of these values are shifted by the largest elevation
    to place the sources and receivers within the model.

    :param model: Layered model object
    :type model: :class:`~pyrocko.cake.LayeredModel
    :param origin: Event origin object
    :type origin: :class:`~obspy.core.event.Origin`
    :param inventory: Station inventory object
    :type inventory: :class:`~obspy.core.inventory.Inventory`
    :param phases: List of phases to model, defaults to ['P','S']
    :type phases: list, optional
    :returns:
     - **results** (*dict*) -- dictionary keyed by NET.STA codes from
        input **inventory** with values that are structured as follows:
        [:meth:`~pyrocko.cake.LayeredModel.arrivals` raw output, e_offset]

         - **e_offset** is the vertical offset applied to modeled arrivals to place
           receivers (stations) and the origin inside the provided velocity model
           This value should be **added** to all elevations (y coordinates) of the
           ray-path in the arrivals raw output to restore the true elevations 
    """    
    results = {}
    # Get Origin Hypocentral parameters
    olat = origin.latitude
    olon = origin.longitude
    odep = origin.depth
    oele = -1.*origin.depth
    time = origin.time

    results = {}
    # Iterate across stations
    for net in inventory.networks:
        for sta in net.stations:
            rlat = sta.latitude
            rlon = sta.longitude
            # Get epicentral distance
            delta = locations2degrees(olat, olon, rlat, rlon)
            # Handle positive elevations
            rele = sta.elevation
            rdep = -1.*rele
            # Get maximum elevation
            max_ele = max(oele, rele)
            # Determine if offset is needed to adjust model datum
            if max_ele > 0:
                d_offset = -1.*max_ele
            else:
                d_offset = 0
            # Calculate arrival times for specified phases
            arrivals = model.arrivals(
                distances=[delta],
                zstart = odep - d_offset,
                zstop = rdep - d_offset,
                phases=phases
            )
            routing = 's2r'
            # If the source-to-receiver ray does not populate
            # try reversing the propagation direction
            if arrivals == []:
                arrivals = model.arrivals(
                    distances=[delta],
                    zstart=rdep - d_offset,
                    zstop=odep - d_offset,
                    phases=phases
                )
                routing = 'r2s'
            results.update({f'{net.code}.{sta.code}': [arrivals, d_offset, routing]})
    return results


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

def ray2pick(ray, wfid, origin, method_str):
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
    mid = method_str
    pick = Pick(time=ta,
                phase_hint=ray.used_phase().given_name(),
                waveform_id=wfid,
                evaluation_mode='automatic',
                method_id=ResourceIdentifier(id=mid),
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
                    _arr.used_phase().given_name(),
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