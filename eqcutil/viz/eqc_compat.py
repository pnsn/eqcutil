"""
:module: eqcorrscan_utils.snuffler
:auth: Nathan T. Stevens
:email: ntsteven@uw.edu
:org: Pacific Northwest Seismic Network
:license: GNU GPLv3
:purpose:
    This module extends the pyrocko.obspy_compat.plant method
    to also allow EQcorrscan Template and Tribe objects to be
    visualized with `snuffler`, automaticaly including the
    event information as marker(s) and scaling the maximum
    number of visible traces to the number of unique trace ID's
    present in either object's waveform data.
"""
from collections import defaultdict
from obspy import Catalog, Stream, UTCDateTime, read_inventory, Inventory
from obspy.core.event import *
from pyrocko import obspy_compat, model
from pyrocko.gui.snuffler.marker import Marker, EventMarker, PhaseMarker

def plant():
    """
    Runs the :meth:`~pyrocko.obspy_compat.plant` method and extends
    it to also add the `snuffle` method to EQcorrscan Template and Tribe
    class objects.
    """    
    obspy_compat.plant()
    import eqcorrscan
    eqcorrscan.Tribe.snuffle = snuffle_tribe
    eqcorrscan.Template.snuffle = snuffle_template

def pick_to_phase(pick, hash=None, kind=0):
    """Convert an obspy Pick object into a Snuffler PhaseMarker object

    :param pick: pick object
    :type pick: obspy.core.event.Pick
    :param hash: event hash to associate with this pick, defaults to None
        This comes from Snuffler Event objects
    :type hash: str, optional
    :param kind: integer kind code (color) to assign this marker, defaults to 0
        Accepted values: 0-6
    :type kind: int, optional
    :return: phase marker
    :rtype: ~.PhaseMarker
    """    
    if pick.evaluation_mode == 'automatic':
        automatic=True
    else:
        automatic=False
    tp = pick.time
    dt = pick.time_errors['uncertainty']
    if isinstance(dt, float):
        tmin = tp - dt
        tmax = tp + dt
    else:
        tmin = tp
        tmax = tp
    
    if hasattr(pick, 'phase_hint'):
        phase_hint = pick.phase_hint
    else:
        phase_hint=None

    pmarker = PhaseMarker(
        tmin=tmin.timestamp,
        tmax=tmax.timestamp,
        nslc_ids=[tuple(pick.waveform_id.id.split('.'))],
        kind=kind,
        event_hash=hash,
        phasename=phase_hint,
        automatic=automatic
    )
    return pmarker

def phase_to_pick(phase):
    """Convert a Snuffler PhaseMarker to an obspy Pick object

    If the phase marker stretches across a range of time
    i.e., tmin < tmax, then the pick time is set as the
    window-centered time and a time-error is attached to the
    output pick object scaled to the half-width of the window

    :param phase: phase marker
    :type phase: ~.PhaseMarker
    :return: pikc object
    :rtype: obspy.core.event.Pick
    """    
    if hasattr(phase,'automatic'):
        evaluation_mode = 'automatic'
    else:
        evaluation_mode = 'manual'
    
    tmin = phase.tmin
    tmax = phase.tmax
    if tmin == tmax:
        dt = None
        tp = UTCDateTime(tmin)
    else:
        dt = 0.5*(tmax - tmin)
        tp = UTCDateTime(tmin) + dt

    nslc = '.'.join(list(phase.nslc_ids[0]))
    
    pick = Pick(
        resource_id=ResourceIdentifier(prefix='smi:local/eqc_compat/phase_to_pick'),
        time=tp,
        time_errors=QuantityError(uncertainty=dt),
        waveform_id=WaveformStreamID(seed_string=nslc),
        evaluation_mode = evaluation_mode,
        phase_hint=phase.get_phasename()
        )
    return pick

def to_pyrocko_events_and_picks(catalog, altnames=None):
    """Convert an obspy Catalog object with picks and origins into 
    lists of Snuffler Event objects, and EventMarkers and PhaseMarker objects
    with the option to assign alternative names to each event.

    :param catalog: event catalog with origins and picks
        Should be only one origin per event, or events should have
        their preferred_origin_id attribute set
    :type catalog: obspy.core.event.Catalog
    :param altnames: list of alternative names for each event in `catalog`, defaults to None
    :type altnames: list-like, optional
    :returns: 
    - **events** (*list* of *~.Event*) - list of Snuffler Event objects
    - **markers** (*list* of *~.PhaseMarker* and *~.EventMarker*) - list of Snuffler Event and Phase markers
        that are associated
    """    
    ocat = catalog
    if ocat is None:
        return None
    
    events = []
    markers = []
    for _e, oevent in enumerate(ocat):
        orig = oevent.preferred_origin()
        if orig is None:
            orig = oevent.origins[0]
       
        if altnames is None:
            name = f'{oevent.resource_id}'
        else:
            name = altnames[_e]

        event = model.Event(name=name,
                            time=orig.time.timestamp,
                            lat=orig.latitude,
                            lon=orig.longitude,
                            depth=orig.depth,
                            region=orig.region)
        events.append(event)
        emarker = EventMarker(event=event)
        markers.append(emarker)
        hash = emarker.get_event_hash()
        for pick in oevent.picks:
            phase = pick_to_phase(pick, hash=hash)
            markers.append(phase)
    return events, markers

def snuffle_picked_catalog(st, catalog, altnames=None, inventory=None):
    events, markers = to_pyrocko_events_and_picks(catalog, altnames=altnames)
    if isinstance(inventory, str):
        inv = read_inventory(inventory)
    elif isinstance(inventory, Inventory):
        inv = inventory
    else:
        inv = None
    
    nslc_set = set([tr.id for tr in st])

    exit_code, marker_pile = st.snuffle(markers=markers, ntracks=len(nslc_set), inventory=inv)
    return exit_code, marker_pile


def markers_to_cat(markers):
    cat = Catalog()
    mdict = defaultdict(list)
    for m in markers:
        _hash = m.get_event_hash()
        if _hash is None:
            _hash = 'unassoc'
        mdict[_hash].append(m)
    
    for _k, _v in mdict.items():
        pmark = []
        for _m in _v:
            if isinstance(_m, EventMarker):
                emark = _m
                ev= _m.get_event()
                origin = Origin(resource_id=ResourceIdentifier(id=f'smi:local/{_m.get_event_hash()}'),
                                time=UTCDateTime(ev.time),
                                latitude=ev.lat,
                                longitude=ev.lon,
                                depth=ev.depth)
                if ev.name == 'Event':
                    event = Event(resource_id=ResourceIdentifier(prefix='smi:local/new_event'))
                else:
                    event = Event(resource_id=ResourceIdentifier(id=ev.name))
                event.origins.append(origin)
                event.preferred_origin_id = origin.resource_id
                
            elif isinstance(_m, PhaseMarker):
                pmark.append(_m)
        
        if _k == 'unassoc':
            event = Event(resource_id=ResourceIdentifier(prefix='smi:local/unassociated'))

        for p in pmark:
            pick = phase_to_pick(p)
            if _k != 'unassoc':
                assoc = Arrival(pick_id=pick.resource_id, phase=p.get_phasename())
                event.preferred_origin().arrivals.append(assoc)
            event.picks.append(pick)

        
        cat.events.append(event)
    return cat

        

        




def snuffle_template(template, altname=None, **kwargs):
    """
    Initialize a Pyrocko "snuffler" GUI minimally displaying the
    traces and event marker contained in this :class:`~eqcorrscan.core.match_filter.template.Template`
    object's **st** and **event** attributes.

    Note: adding an inventory for relevant stations will allow distance sorting!

    :return: (return_tag, markers)
    :rtype: (str, list[pyrocko markers] )
    """
    if 'ntracks' not in kwargs.keys():
        kwargs.update({'ntracks': len({tr.id for tr in template.st})})
    cat = Catalog()
    cat.events.append(template.event)
    # If a supplementary catalog is provided
    if 'catalog' in kwargs.keys():
        # If it is not identical to the catalog comprising all templates' events
        if kwargs['catalog'] != cat:
            cat += kwargs.pop('catalog')
        # Skip if identical (pop to nowhere)
        else:
            kwargs.pop('catalog')
    if altname is None:
        altname = template.name
    else:
        altname = str(altname)
    events, markers = to_pyrocko_events_and_picks(cat, altname=altname)
    kwargs.update({'markers': markers})
    return template.st.snuffle(**kwargs)

def snuffle_tribe(tribe, altnames=None, **kwargs):
    """Initialize a Pyrocko "snuffler" GUI instance
    for the contents of this :class:`~eqcorrscan.core.match_filter.tribe.Tribe`
    minimally displaying waveform data for all **template.st** and event
    markers for all **template.event** in this tribe

    :param tribe: tribe to snuffle
    :type tribe: eqcorrscan.core.match_filter.tribe.Tribe
    :return: (return_tag, markers)
    :rtype: (str, list[pyrocko markers] )
    """    
    big_st = Stream()
    big_markers = []

    # Get events and associated event markers
    if 'preferred' in kwargs.keys():
        pref = kwargs.pop('preferred')
    else:
        pref = True

    if altnames is None:
        altnames = [tmp.name for tmp in tribe]

    for _e, temp in enumerate(tribe):
        big_st += temp.st.copy()
        events, markers = to_pyrocko_events_and_picks(
            Catalog(events=[temp.event.copy()]),
            altname=altnames[_e],
            preferred=pref)
        big_markers += markers    
    
    if 'markers' in kwargs.keys():
        if all(isinstance(_e, Marker) for _e in kwargs['markers']):
            big_markers += kwargs['markers']

    kwargs.update({'markers': big_markers})
    # # If a supplementary catalog is provided
    # if 'catalog' in kwargs.keys():
    #     # If it is not identical to the catalog comprising all templates' events
    #     if kwargs['catalog'] != cat:
    #         cat += kwargs.pop('catalog')
    #     # Skip if identical (pop to nowhere)
    #     else:
    #         kwargs.pop('catalog')

    # If ntracks is not provided, set ntracks to the number of unique channel IDs
    if 'ntracks' not in kwargs.keys():
        kwargs.update({'ntracks': len({tr.id for tr in big_st})})

    return big_st.snuffle(**kwargs)


