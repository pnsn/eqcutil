import numpy as np
import pandas as pd

from eqcorrscan import Template

def _rms(array):
    """
    Calculate RMS of array.

    :type array: numpy.ndarray
    :param array: Array to calculate the RMS for.

    :returns: RMS of array
    :rtype: float
    """
    return np.sqrt(np.mean(np.square(array)))


def get_template_pick_snrs(template, prepick_scalar=0.9):
    if not isinstance(template, Template):
        raise TypeError('template must be type eqcorrscan.Template')
    if not isinstance(prepick_scalar, float):
        raise TypeError('prepick_scalar must be type float')
    elif not 0 < prepick_scalar <= 1.:
        raise ValueError('prepick scalar must be in range (0, 1]')
    
    snr_out = []
    for tr in template.st:
        id = tr.id
        full_rms = _rms(tr.copy().data)
        noise_tr = tr.copy().trim(endtime=tr.stats.starttime + template.prepick*prepick_scalar)
        signal_tr = tr.copy().trim(starttime=tr.stats.starttime + template.prepick,
                                   endtime = tr.stats.starttime + (1. + prepick_scalar)*template.prepick)
        noise_rms = _rms(noise_tr.data)
        signal_rms = _rms(signal_tr.data)
        snr_out.append([id, full_rms, noise_rms, signal_rms, signal_rms/noise_rms])
    out = pd.DataFrame(snr_out, columns=['trace_id','full_trace_rms','noise_rms','signal_rms','snr'])
    return out
