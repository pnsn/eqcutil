"""
:module: eqcorrscan_utils.catalog
:auth: Nathan T. Stevens; Barret Johnson
:email: ntsteven@uw.edu; bnjo@wu.edu
:org: Pacific Northwest Seismic Network
:license: GNU GPLv3
:purpose: This module contains methods that help with
    metadata manipulation in :class:`~obspy.core.event.Catalog`
    objects
"""
import logging
import obspy
import numpy as np
from obspy.core.event import WaveformStreamID

import eqcorrscan.utils.catalog_utils as eqcu
from eqcorrscan.utils.clustering import catalog_cluster

Logger = logging.getLogger(__name__)

def apply_phase_hints(catalog):
    """Apply phase hints from arrivals to picks in this catalog
    if they are used in the preferred origin of each event in
    this catalog.

    Note: Changes are made in-place, altering the input. If you
    want to preserve your data, use :meth:`~obspy.core.event.Catalog.copy()`
    on your input.

    :param catalog: input catalog
    :type catalog: :class:`~obspy.core.event.Catalog
    :return:
     - **catalog** (*obspy.core.event.Catalog*) -- updated catalog
    """    
    if not isinstance(catalog, obspy.core.event.Catalog):
        raise TypeError('catalog must be type obspy.core.event.Catalog')
    for event in catalog.events:
        # try: 
        origin = event.preferred_origin()
        # except ????: Make sure this isnt hanging
        #     event = event.origins[0]
        arrivals = origin.arrivals
        # Iterate across arrivals and ensure
        for _arr in arrivals:
            # Get 
            phz = _arr.phase
            # Get pick
            pick = _arr.pick_id.get_referred_object()
            # If the pick is present
            if pick != None:
                pick.phase_hint = str(phz)
            else:
                msg = f'Missing pick for ORID: {origin.resource_id} & ARID: {_arr.resource_id}'
                Logger.warning(msg)
    return catalog


def filter_picks(catalog, **kwargs):
    """Extension of :meth:`~eqcorrscan.utils.catalog_utils.filter_picks`, adding
    an option 'preferred' for argument **enforce_single_pick**. This option
    Uses the arrivals from the preferred origin of each event in **catalog**
    to filter out extraneous picks and then applies filtering specified as normal

    Filter events in the catalog based on a number of parameters.

    :param catalog: Catalog to filter.
    :type catalog: obspy.core.event.Catalog
    :param stations: List for stations to keep picks from. Defaults to None
    :type stations: list or Nonetype, optional
    :param channels: List of channels to keep picks from. Defaults to None
    :type channels: list or NoneType, optional
    :param networks: List of networks to keep picks from. Defaults to None
    :type networks: list or NoneType, optional
    :param locations: List of location codes to use. Defaults to None
    :type locations: list or NoneType, optional
    :param top_n_picks: Filter only the top N most used station-channel pairs.
        Defaults to None
    :type top_n_picks: int or NoneType, optional
    :param evaluation_mode:
        To select only 'manual' or 'automatic' picks, or use 'all'. Defaults to 'all'.
    :type evaluation_mode: str, optional
    :param phase_hints: List of retained phase hints, or None to use all. Defaults to None
    :type phase_hints: list or NoneType, optional
    :param enforce_single_pick:
        Method to enforce using only one pick of each phase-hint per
        station or False to leave all. Can be {False, "earliest", "preferred"}.
        Defaults to False
    :type enforce_single_pick: str or bool

    :return:
        Filtered Catalog - if events are left with no picks, they are removed
        from the catalog.
    :rtype: obspy.core.event.Catalog

    .. note::
        Will filter first by station, then by channel, then by network, if
        using top_n_picks, this will be done last, after the other filters
        have been applied.

    .. note::
        Doesn't work in place on the catalog, your input catalog will be safe
        unless you overwrite it.

    .. note:: Doesn't expand wildcard characters.
    """    
    if 'enforce_single_pick' in kwargs.keys():
        if kwargs['enforce_single_pick'] in [False,'earliest']:
            cat = eqcu.filter_picks(catalog, **kwargs)
        elif kwargs['enforce_single_pick'] == 'preferred':
            kwargs.update({'enforce_single_pick': False})
            for event in catalog:
                origin = event.preferred_origin()
                ppicks = [_arr.pick_id.get_referred_object() for _arr in origin.arrivals]
                event.picks = ppicks  
            # Force only preferred
            event.origins = [event.preferred_origin()]
            cat = eqcu.filter_picks(catalog, **kwargs)

    else:
        cat = eqcu.filter_picks(catalog, **kwargs)
    return cat


def transfer_picks(catalog, mapping={}):
    if len(mapping) > 0:
        for _e, event in enumerate(catalog.events):
            for _p, pick in enumerate(event.picks):
                seed = pick.waveform_id.id
                if seed in mapping.keys():
                    pick.waveform_id = WaveformStreamID(seed_string=mapping[seed])
    return catalog

# def filter_picks(catalog, stations=None, channels=None, networks=None,
#                  locations=None, top_n_picks=None, evaluation_mode='all',
#                  phase_hints=None, enforce_single_pick=False,
#                  min_delta=None, max_delta=None, min_baz=None, max_baz=None,
#                  min_depth=None, max_depth=None):
#     """Filter picks in an ObsPy :class:`~obspy.core.event.Catalog` object based
#     on station, phase-type, and source-receiver geometries.
    
#     Augments the EQcorrscan :meth:`~eqcorrscan.utils.catalog_utils.filter_picks` 
#     method by adding the **min_delta**, **max_delta**, **min_az**, **max_az**,
#     **min_depth**, **max_depth**, and **inv** arguments to facilitate
#     source-receiver geometry based filtering of :class:`~obspy.core.event.Pick`
#     objects in the input :class:`~obspy.core.event.Catalog` object for picks
#     that are associated to one or more :class:`~obspy.core.event.Origin` objects
#     via :class:`~obspy.core.event.Arrival` objects.

#     Parameters
#     ----------
#     :param catalog: Catalog to filter.
#     :type catalog: obspy.core.event.Catalog
#     :param stations: List for stations to keep picks from.
#     :type stations: list
#     :param channels: List of channels to keep picks from.
#     :type channels: list
#     :param networks: List of networks to keep picks from.
#     :type networks: list
#     :param locations: List of location codes to use
#     :type locations: list
#     :param top_n_picks: Filter only the top N most used station-channel pairs.
#     :type top_n_picks: int
#     :param evaluation_mode:
#         To select only manual or automatic picks, or use all (default).
#     :type evaluation_mode: str
#     :param phase_hints: List of retained phase hints, or None to use all
#     :type phase_hints: list
#     :param enforce_single_pick:
#         Method to enforce using only one pick of each phase-hint per
#         station or False to leave all. Can be {False, "earliest"}
#     :type enforce_single_pick: str
#     :param min_delta: minmum great-circle distance for source-receiver separation
#         in degrees, defaults to 0.
#     :type min_delta: float, optional.
#     :param max_delta: maximum great-circle distances for source-receiver separation
#         in degrees, defaults to 180.
#     :return:
#      - **out** (*obspy.core.event.Catalog*) -- filtered catalog
#     """
#     spatial_kwargs = {'min_delta': min_delta,
#                       'max_delta': max_delta,
#                       'min_baz': min_baz,
#                       'max_baz': max_baz,
#                       'min_depth': min_depth,
#                       'max_depth': max_depth}
#     # Handle Nones
#     for _k, _v in spatial_kwargs:
#         if _v is None:
#             if 'max' in _k:
#                 spatial_kwargs.update({_k: np.inf})
#             else:
#                 spatial_kwargs.update({_k: -np.inf})
    
#     # Run filtering
    

#     # Run distance filtering first
#     if min_delta > 0 or max_delta < 180:
#         for event in catalog:
#             picks = []
#             origin = event.preferred_origin()
#             if len(event.picks) == 0:
#                 continue
#             for arrival in origin.arrivals:
#                 if min_delta <= arrival.distance <= max_delta:
#                     pick = arrival.waveform_id.get_referred_object()
#                     picks.append(pick)
#             event.picks = picks

#     # Then run the rest of the filtering
#     catalog = filter_picks(
#         catalog,stations=stations, channels=channels, networks=networks,
#         locations=locations, top_n_picks=top_n_picks,
#         evaluation_mode=evaluation_mode, phase_hints=phase_hints,
#         enforce_single_pick=enforce_single_pick)

#     return catalog

