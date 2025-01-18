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

from obspy import Catalog, Stream, UTCDateTime
from obspy.core.event import Pick, QuantityError, WaveformStreamID, ResourceIdentifier
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
    if phase.automatic:
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

def to_pyrocko_events_and_picks(catalog, altname=None, preferred=True):
    ocat = catalog
    if ocat is None:
        return None
    
    events = []
    markers = []
    for oevent in ocat:
        if preferred:
            origs = [oevent.preferred_origin()]
        else:
            origs = oevent.origins
        for orig in origs:
            if altname is None:
                name = f'{oevent.resource_id}-{orig.resource_id}'
            else:
                name = altname

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


