import copy, logging, os
from pathlib import Path

Logger = logging.getLogger(__name__)

def rename_templates(tribe, include_contributor=True, inplace=True):
    """Rename templates in a tribe with their PNSN EVID's
    with the option to enact the changes on a deepcopy of the
    input :class:`~eqcorrscan.core.match_filter.tribe.Tribe.

    The name based on the **event.resource_id** attribute
    of each :class:`~eqcorrscan.core.match_filter.template.Template`
    object within the input **tribe**, which should follow
    QuakeML formatting conventions.

    I.e., the last 2 "/" delimited elements are:
    [-2] contributing network code
    [-1] event ID for that contributing network

    E.g., For resource_id = 'quakeml:uw.anss.org/Event/UW/61275716'
    would yield 'uw61275716' as the re-assigned name
    if include_contributor = True. 

    Parameters
    ----------
    :param tribe: tribe containing templates to relabel
    :type tribe: :class:`~eqcorrscan.core.match_filer.Tribe
    :param include_contributor: include contributor ID code
        in lower-case letters? Defaults to True
    :type include_contributor: bool, optional
    :param inplace: should changes be made in-place on 
        the input **tribe**? Defaults to True
    :return:
     - **tribe** (*eqcorrscan.core.match_filter.Tribe*) - tribe of renamed templates
    """
    # If not inplace, make a deepcopy of the tribe within the
    # scope of this method
    if not inplace:
        tribe = copy.deepcopy(tribe)
    for template in tribe:
        event = template.event
        oldname = template.name
        parts = str(event.resource_id).split(sep='/')
        evid = parts[-1]
        contributor = parts[-2]
        if include_contributor:
            name = f'{contributor.lower()}{evid}'
        else:
            name = str(evid)
        Logger.debug(f'renamed template {oldname} -> {name}')
        template.name = name
    return tribe


def augment_template(template, client, padding= 120., min_ncomponents=3):
    """Retrieve additional waveform data for missing channels  for each 
    instrument in a EQcorrscan :class:`~eqcorrscan.core.match_filter.template.Template`
    object's **st** attribute from the provided **client**. If there are
    new trace ID's (NSLC codes) retrieved, they are pre-processed to match
    the sampling rate, bandpass filtering, and timing of the first matching
    trace for the relevant instrument code (NSLC code minus the component character)

    .. rubric:: Explainer
        Say a template has one entry for 'UW.SHUK..BHZ' and the following
        pre-processing steps:
        - 50 Hz sampling rate
        - 0.5 - 20 Hz bandpass filtering (4th order)

        this method will fetch data from the client for 'UW.SHUK..BH?' for 
        the start and end time of the 'UW.SHUK..BHZ' trace--padded by **padding**
        seconds--and apply pre-processing in the following order
        - resample 
            - downsampling uses :meth:`~obspy.core.trace.Trace.resample` with no_filter=False,
            - upsampling uses :meth:`~obspy.core.trace.Trace.interpolate` with method='lanczos')
        - filter
        - trim

    Parameters
    ----------
    :param template: template to augment, warning: modifications are made in-place
    :type template: eqcorrscan.core.match_filter.template.Template
    :param client: client object with a :meth:`~obspy.clients.fdsn.Client.get_waveforms`-type method
    :type client: obspy.clients.fdsn.Client or similar (e.g., obsplus.bank.wavebank.WaveBank)
    :param padding: amount of padding in seconds to add to each retrieved trace
        to capture and exclude filtering edge-effects, defaults to 120.
    :type padding: float-like, optional
    :param min_ncomponents: minimum number of expected components for 
        instruments, defaults to 3
    :type min_ncomponents: int, optional
    :return: augmented template
    :rtype: eqcorrscan.core.match_filter.template.Template
    """    
    # Get unique instrument codes
    insts = {tr.id[:-1] for tr in template.st}
    # Iterate across codes and see if there are missing components
    for inst in insts:
        sub_st = template.st.select(id=inst + '?')
        t0 = sub_st[0].stats.starttime
        t1 = sub_st[0].stats.endtime
        if len(sub_st) < min_ncomponents:
            n, s, l, c = inst.split('.')
            # Fetch waveforms
            ist = client.get_waveforms(network=n, station=s, location=l, channel=c + '?',
                                       starttime=t0 - padding, endtime=t1 + padding)
            # If there are new waveforms
            if len(ist) > len(sub_st):
                # Iterate across all fetched waveforms
                for tr in ist:
                    # If the trace ID is not in sub_st IDs
                    if tr.id not in [tr.id for tr in sub_st]:
                        # Confirm that S/R matches
                        if tr.stats.sampling_rate > template.samp_rate:

                            Logger.info(f'downsampling new trace to match template: {tr.id} {tr.stats.sampling_rate} -> {template.samp_rate}')
                            tr.resample(template.samp_rate, no_filter=False)
                        elif tr.stats.sampling_rate < template.samp_rate:
                            Logger.info(f'upsampling new trace to match template: {tr.id} {tr.stats.sampling_rate} -> {template.samp_rate}')
                            tr.interpolate(method='laczos')
                        else:
                            pass
                        Logger.info(f'Bandpass filtering new trace to match template preprocessing')
                        tr.filter('bandpass',
                                  freqmin = template.lowcut,
                                  freqmax=template.highcut,
                                  corners=template.filt_order)
                        tr.trim(starttime=t0, endtime=t1)
                        template.st += tr
    return template

