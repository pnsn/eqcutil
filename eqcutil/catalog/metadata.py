
import logging
from obspy import Catalog

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
    if not isinstance(catalog, Catalog):
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