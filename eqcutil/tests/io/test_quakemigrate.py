import unittest

from eqcutil.io.quakemigrate import *


class TestStream_ID_Formatter(unittest):

    def setUp(self):
        self.phase = 'P'
        self.badphase = 1
        self.mismatch_phase = 'pP'
        self.station = 'MBW2',
        self.chan_mapping={'P': 'HHZ','S':'HHN'}
        self.bad_chan_mapping = {'p': 'HHZ', 'S':'HHN'}
        self.alt_location='--'

    def tearDown(self):
        del self.phase, self.badphase, self.station
        del self.chan_mapping, self.alt_location

    def test_standard(self):
        rid = stream_id_formatter()