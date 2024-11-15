import os, tarfile

from obspy.clients.fdsn import Client
from obspy.core.event import Event
from obspy import Stream, UTCDateTime

from obsplus import WaveBank

from eqcorrscan import Template
from eqcorrscan.utils.pre_processing import _multi_filter
from eqcorrscan.core.match_filter.helpers import _safemembers


class ClientTempate(Template):

    def __init__(self, wave_client, event, name_elements=2, padding=5., prepick=2., length=30., samp_rate=None, lowcut=None, highcut=None, filt_order=None, process_length=None):

        # Ensure that the client provided has this method
        if not hasattr(wave_client, 'get_waveforms_bulk'):
            raise AttributeError(f'wave_client does not have method `get_waveforms`')
        # Handle obsplus wavebank input
        elif isinstance(wave_client, WaveBank):
            self.client = wave_client
            self.client_class = 'obsplus.bank.wavebank.WaveBank'
            self.client_params = {'base_path': wave_client.bank_path,
                                  'path_structure': wave_client.path_structure,
                                  'name_structure': wave_client.name_structure}
        # Handle obspy client input
        elif isinstance(wave_client, Client):
            self.client = wave_client
            self.client_class = 'obspy.clients.fdsn.Client'
            self.client_params = {'base_url': self.client.base_url}

        else:
            raise NotImplementedError('Currently supports obspy.clients.fdsn.Client and obsplus.bank.wavebank.WaveBank')
        
        if not isinstance(event, Event):
            raise TypeError('event must be type obspy.core.event.Event')


        if not isinstance(padding, (int, float)):
            raise TypeError
        elif padding < 0:
            raise ValueError
        else:
            self.padding = float(padding)
        
        if not isinstance(prepick, (int, float)):
            raise TypeError
        elif prepick < 0:
            raise ValueError
        else:
            prepick = float(prepick)

        if not isinstance(length, (int, float)):
            raise TypeError
        elif length <= 0:
            raise ValueError
        else:
            self.length = float(length)

        self.bulk = self._compose_bulk_from_picks()
        st = self._st_from_bulk(self, self.bulk)

        super().__init__(self, event=event, st=st, name=''.join(event.resource_id.split('/')[-1*name_elements:]), 
                         lowcut=lowcut, highcut=highcut, samp_rate=samp_rate, filt_order=filt_order,
                         prepick=prepick, process_length=process_length)
        
        


    def _compose_bulk_from_picks(self):
        bulk = []
        for pick in self.event.picks:
            line = pick.waveform_id.id.split('.')
            line += [pick.time - self.prepick - self.padding, pick.time + self.length + self.padding]
            line = tuple(line)
            bulk.append(line)
        return bulk
    
    def _write_bulk(self, dirname):
        with open(os.path.join(dirname, 'bulk.csv'), 'w') as file:
            for line in self.bulk:
                for _e, _p in enumerate(line):
                    file.write(f'{_p}')
                    if _e < 5:
                        file.write(',')
                    else:
                        file.write('\n')

    def _read_bulk(self, filename):
        path, file = os.path.split(filename)
        if file != 'bulk.csv':
            raise ValueError
        else:
            bulk = []
            with open(filename, 'r') as file:
                lines = file.readlines()
            for line in lines:
                elements = line.split(',')
                elements[4] = UTCDateTime(elements[4])
                elements[5] = UTCDateTime(elements[5])
                bulk.append(tuple(elements))
        return bulk
        
    def _st_from_bulk(self, bulk, max_workers=4):
        # Make sure we don't exceed the needed/available resources
        max_workers = min(len(st), max_workers, os.cpu_count())
        # Get waveforms from client
        st = self.client.get_waveforms_bulk(bulk)
        # Apply filtering
        st = _multi_filter(st, self.highcut, self.lowcut, self.filt_order,
                           max_workers=max_workers, chunksize = len(st)//max_workers)
        # Trim off padding
        for tr in st:
            tr.trim(starttime=tr.stats.startime + self.padding,
                    endtime=tr.stats.endtime - self.padding)
        return st


    def write(self, dirname, compress=True):
        # Construct template-specific directory
        dirname = os.path.join(dirname, self.name)
        if not os.path.isdir(dirname):
            os.makedirs(dirname)

        # Write client parameters
        with open(os.path.join(dirname, 'client_parameters.csv'), 'w') as file:
            for _k, _v in self.client_params.items():
                file.write(f'{key}: {str(_v)}, ')

        # Write most other parameters
        with open(os.path.join(dirname, 'template_parameters.csv'), 'w') as file:
            for key in self.__dict__.keys():
                if key not in ['st','event','client','client_params','name', 'st_bulk']:
                    file.write(f'{key}: {str(self.__dict__[key])}, ')
        # Write Event
        self.event.write(os.path.join(dirname, 'event.xml'), format='QuakeML')
        # Write Bulk
        self._write_bulk(dirname)

        if compress:
            if not dirname.endswith('.tgz'):
                dirname += '.tgz'
            with tarfile.open(dirname, "w:gz") as tar:
                tar.add(dirname, arcname=os.path.basename(dirname))

    def _read_from_folder(self, dirname):
        

    def read(self, filename):
        fname, ext = os.path.splitext(filename)
        if ext == '.tgz':
            with tarfile.open(filename, "r:*") as arc:
                temp_dir = tempfile.mkdtemp()
                arc.extractall(path=temp_dir, members=_safemembers(arc), filter=None)
        elif ext == '':
            self._read_from_folder(dirname=filename)

            