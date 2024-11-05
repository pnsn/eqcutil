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

from obspy import Catalog, Stream, Inventory
from pyrocko import obspy_compat

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

def snuffle_template(template, **kwargs):
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
    if 'catalog' not in kwargs.keys():
        kwargs.update({'catalog': Catalog(events=[template.event])})
    return template.st.snuffle(catalog=Catalog([template.event]), **kwargs)

def snuffle_tribe(tribe, **kwargs):
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
    cat = Catalog()
    for temp in tribe:
        big_st += temp.st
        cat.events.append(temp.event)

    # If ntracks is not provided, set ntracks to the number of unique channel IDs
    if 'ntracks' not in kwargs.keys():
        kwargs.update({'ntracks': len({tr.id for tr in big_st})})
    
    # If a catalog is provided
    if 'catalog' in kwargs.keys():
        # If it does not perfectly match the catalog comprised of template events
        if kwargs['catalog'] != cat:
            # Add them together
            kwargs['catalog'] += cat
        # Otherwise, catalog provided matches, don't duplicate
        else:
            pass
    # If no catalog is provided, add this to kwargs
    else:
        kwargs.update({'catalog': cat})
    # Run instance and capture output
    return big_st.snuffle( **kwargs)

