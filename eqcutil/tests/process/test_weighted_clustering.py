from unittest import TestCase

import numpy as np

from obspy.clients.fdsn import Client
from obspy import read_events, UTCDateTime, Stream

from eqcutil.process.weighted_clustering import (
    _correlation_preprocessing,
    _run_normxcorr,
    _get_coherances,
    _remap_coherances_to_secondaries
)

def fetch_example_templates():
    client = Client('IRIS')
    event1 = {'evid':'uw10738413',
              'etype': 'lf',
              't0': UTCDateTime('2008-03-15T04:28:20.42'),
              'stas':['UW.MBW.*.EHZ','UW.RPW.*.EHZ']}
    event2 = {'evid': 'uw10772278',
              'etype': 'eq',
              't0': UTCDateTime('2009-07-14T04:31:07.59'),
              'stas': ['UW.MBW.*.EHZ','UW.PASS..BHZ']}
    event3 = {'evid': 'uw61280137',
              'etype': 'lf',
              't0': UTCDateTime('2017-07-07T13:19:32.15'),
              'stas': ['UW.MWB.*.EHZ','UW.RPW.*.EHZ','UW.SHUK..BHZ']}
    event4 = {'evid': 'uw10736828',
              'etype':'eq',
              't0': UTCDateTime('2007-09-29T19:42:34.65'),
              'stas': ['UW.MBW.*.EHZ','UW.RPW.*.EHZ']}
    for _e in [event1, event2, event3, event4]:
        bulk = []
        for _sta in _e['stas']:
            line = tuple(_sta.split('.') + [_e['t0'] - 5, _e['t0']+ 50])
            bulk.append(line)
        st = client.get_waveforms_bulk(bulk)
        _e.update({'st': st})
    return event1, event2, event3, event4




class TestWeigtedClustering(TestCase):
    e1,e2,e3, e4 = fetch_example_templates()
    def setUp(self):
        self.primary = self.e1['st'].copy().filter('bandpass', freqmin=0.2, freqmax=20).resample(50)
        self.secondaries = [_e['st'].copy().filter('bandpass', freqmin=0.2, freqmax=20).resample(50) for _e in [self.e2, self.e3, self.e4]]

    def tearDown(self):
        del self.primary
        del self.secondaries

    def test__correlation_preprocessing(self):
        sbackup = [_st.copy() for _st in self.secondaries]
        ppp, pps, idc = _correlation_preprocessing(self.primary.copy(), sbackup, shift_len=5.0)
        self.assertIsInstance(ppp, Stream)
        self.assertIsInstance(pps, list)
        for _e in pps:
            self.assertIsInstance(_e, Stream)
        self.assertEqual(len(idc), len(pps))
        # Assert that non-matching seed ids have all-nan valued data vectors (but not masked)
        for _e, tmp in enumerate(pps):
            _s = self.secondaries[_e]
            _sid = {_tr.id for _tr in _s}
            for _tr in tmp:
                if _tr.id not in _sid:
                    self.assertTrue(all(~np.isfinite(_tr.data)))
                    self.assertFalse(np.ma.is_masked(_tr.data))
    
    def test__run_normxcorr(self):
        ppp, pps, idcs = _correlation_preprocessing(self.primary, self.secondaries, shift_len=0.0)
        outs = _run_normxcorr(pps, ppp)
        breakpoint()
        


        


