import logging
from obspy import Stream
from eqcorrscan.utils.clustering import *
from eqcorrscan.utils.pre_processing import _prep_data_for_correlation

Logger = logging.getLogger(__name__)



def _correlation_preprocessing(primary: Stream, secondaries: list, shift_len=0.0) -> tuple[Stream, list, list]:
    sr_1 = primary[0].stats.sampling_rate
    n_2 = len(secondaries)
    end_trim = int((shift_len * sr_1)/2)
    # Trim off extra samples that would slide outside the bounds of the primary
    _streams = []
    if end_trim > 0:
        for st in secondaries:
            _st = st.copy()
            for _tr in _st:
                _tr.data = _tr.data[end_trim: -end_trim]
                if _tr.stats.sampling_rate != sr_1:
                    raise NotImplementedError('Sampling rates differ')
            _streams.append(_st)
        secondaries = _streams
    # If no shifts allowed, just make a copy of the list
    else:
        secondaries = [st.copy() for st in secondaries]
    
    # Run prep data for correlation
    primary_pp, secondaries_pp, indices = _prep_data_for_correlation(
        stream=primary.copy(), templates=secondaries,
        template_names = list(range(len(secondaries))), force_stream_epoch=False)
    
    return (primary_pp, secondaries_pp, indices)

def _run_normxcorr(primary_pp: Stream, secondaries_pp: list,
                   xcorr_func='fftw', concurrency='concurrent',
                   cores=1, stack=False, **kwargs) -> tuple[np.ndarray, np.ndarray, list]:
    """Run a normalized multi-channel cross-correlation between a primary template
    and secondary templates with a specified function

    :param primary_pp: primary / reference template
    :type primary_pp: Stream
    :param secondaries_pp: list of secondary / comparison templates
    :type secondaries_pp: list of Stream
    :param xcorr_func: name or method of cross correlation function to use.
        Defaults to 'fftw'
    :type xcorr_func: str or method, optional
    :param concurrency: concurrency mode, defaults to 'concurrent'
    :type concurrency: str, optional
    :param cores: _description_, defaults to 1
    :type cores: int, optional
    :param stack: _description_, defaults to False
    :type stack: bool, optional
    :returns:
     - **cccsums** (*numpy.ndarray*) -- no. secondaries x no. primary channels x no. sample shifts
        array containing normalized cross correlation coefficient values
            Thus, the correlation coefficients for the sliding cross-correlation of the 0th secondary template
            on the 2nd channel specified by the order of primary channels is given as cccsums[0,2,:]
     - **no_chans** (*numpy.ndarray*) -- vector indicating the number of matching traces between
        the input secondary templates and the primary template
     - **matching_chans** (*list*) -- nested list of (STATION, CHANNEL) tuples indicating matching
        station/channel codes between the primary template and secondary templates
    :rtype: tuple[np.ndarray, int, ]
    """    
    mcxc_fun = get_stream_xcorr(xcorr_func, concurrency=concurrency)
    cccsums, no_chans, matching_chans = mcxc_fun(templates=primary_pp, stream=secondaries_pp, cores=cores, stack=stack, **kwargs)
    return (cccsums, no_chans, matching_chans)



def _get_coherances(cccsums, no_chans, samp_rate, shift_len, allow_individual_trace_shifts=False):
    """Calculated the summed normalized correlation coefficients

    :param cccsums: _description_
    :type cccsums: _type_
    :param no_chans: _description_
    :type no_chans: _type_
    :param samp_rate: _description_
    :type samp_rate: _type_
    :param shift_len: _description_
    :type shift_len: _type_
    :param allow_individual_trace_shifts: _description_, defaults to False
    :type allow_individual_trace_shifts: bool, optional
    :return: _description_
    :rtype: _type_
    """    
    sr_1 = samp_rate
    end_trim = int((shift_len * sr_1)/2)
    # If Cross Correlation Coefficients included shifts
    if allow_individual_trace_shifts:
        # Find maximum for each sliding cross correlation, collapsing
        max_ccc_by_channel = cccsums.max(axis=-1) # R3 -> R2
        summed_max_ccc_by_template = max_ccc_by_channel.sum(axis=-1) # R2 -> R1
        coherances = summed_max_ccc_by_template / no_chans # R1 -> R1
        coherances = cccsums.max(axis=-1).sum(axis=-1) / no_chans
    else:
        cccsums = cccsums.sum(axis=1) # R3 -> R2
        coherances = cccsums.max(axis=-1) / no_chans #R2 -> R1
    # Get position in seconds
    positions = (cccsums.argmax(axis=-1) - end_trim) / sr_1

    return (coherances, positions)

def _remap_coherances_to_secondaries(primary_pp, secondaries_pp, indices, coherances, positions, allow_individual_trace_shifts=False):
    _coh = np.empty(len(secondaries_pp))
    if allow_individual_trace_shifts:
        n_max_traces = max([len(_st) for _st in secondaries_pp])
        for _e, _tr in enumerate(primary_pp):
            if np.ma.is_masked(_tr.data):
                positions[:, _e] = np.nan
    # Get the number of shifts from the 1-axis of the positions matrix
    n_shifts_per_stream = positions.shape[1]
    _pos = np.empty([len(secondaries_pp), n_max_traces])
    _pos.fill(np.nan)
    _coh.fill(np.nan)
    # Insert correlations and shifts at the template indices
    _coh[np.ix_(indices)] = coherances
    _pos[np.ix_(indices, range(n_shifts_per_stream))] = (
        positions
    )
    # if not allow_individual_trace_shifts:
    #     _pos = _pos[:, ]



def _get_weighted_coherances(cccsums, matching_chans, primary_pp, chan_weights, shift_len, allow_individual_trace_shifts=False):
    sr_1 = primary_pp.stats.sampling_rate
    end_trim = int((shift_len * sr_1)/2)

    # Confirm that all matching_chans show up in chan_weights
    # Iterate across used station channel tuples
    sc_index = {(_tr.stats.station, _tr.stats.channel): _e for _e, _tr in enumerate(primary_pp)}
    if sc_index.keys() <= chan_weights.keys():
        pass
    else:
        raise KeyError('not all primary station-channel tuples show up in chan_weights keys')
    # Iterate across template pairs
    for _ii, _sc_list in enumerate(matching_chans):
        # Initialize weight sum
        _wgt_sum = 0.
        # Iterate across matched channels
        for _sc in _sc_list:
            # Get the index for appropriate weighting
            _jj = sc_index[_sc]
            # Get the appropriate weight for this channel
            _wgt = chan_weights[_sc]
            # Apply weight to cccsums
            cccsums[_ii, _jj, :] *= _wgt
            # Add _wgt to _wgt_sum
            _wgt_sum += _wgt
    # Safety check for non-positive weight sum
    if _wgt_sum <= 0:
        raise ValueError('Weight summation produced a non-positive weight sum')

    # Proceed as (almost) normal for unweighted coherance
    if allow_individual_trace_shifts:
        coherances = cccsums.max(axis=-1).sum(axis=-1) / _wgt_sum
    else:
        cccsums = cccsums.sum(axis=1)
        coherances = cccsums.max(axis=-1) / _wgt_sum
    positions = (cccsums.argmax(axis=-1) - end_trim) / sr_1
    return (coherances, positions) 

def _remap_coh_and_pos(coherances, positions, primary_pp)





def weighted_cross_chan_correlation(primary, secondaries, shift_len=0.0, channel_weights=1.,
                                    allow_individual_trace_shifts=True, xcorr_func='fftw',
                                    concurrency="concurrent", cores=1, **kwargs):
    minimum_cw_keys = {(_tr.stats.station, _tr.stats.channel) for _tr in primary}
    # For float/int formatted channel_weights
    if isinstance(channel_weights, (int, float)):
        if channel_weights <= 0:
            raise ValueError('float-like channel_weights must be positive-valued')
        cw = {}
        # Populate (s,c) keys from primary
        for _tr in primary:
            cw.update({(_tr.stats.station, _tr.stats.channel): channel_weights})
        
    # For dictionary format channel_weights
    elif isinstance(channel_weights, dict):
        for _k, _v in channel_weights.items():
            if not isinstance(_k, tuple):
                raise TypeError('channel_weights keys must be type tuple')
            elif len(_k) != 2:
                raise ValueError('channel_weights keys must be 2-element tuples')
            elif not all(isinstance(_e, str) for _e in _k):
                raise TypeError('channel_weight key elements must be type str')
            else:
                pass
        cw = channel_weights
    # Check that minimum keys are present
    if minimum_cw_keys <= cw.keys():
            pass
    else:
        raise ValueError('not all (station, channel) code tuples in primary are present in channel_weights')
    # Preprocess templates
    ppp, pps, indces = _correlation_preprocessing(primary, secondaries, shift_len=shift_len)
    # Run cross-correlation
    cccsums, no_matches, matching_chans = _run_normxcorr(ppp, pps, xcorr_func=xcorr_func, concurrency=concurrency,cores=cores, **kwargs)
    # Apply weights & get normalized ccc
    coh, pos = _get_weighted_coherances(cccsums, no_matches, matching_chans, ppp, cw, shift_len, allow_individual_trace_shifts=allow_individual_trace_shifts)
