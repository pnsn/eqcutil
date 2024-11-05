"""
:module: demo.py
:auth: Nathan T. Stevens
:email: ntsteven@uw.edu
:org: Pacific Northwest Seismic Network
:license: GNU GPLv3
:purpose: This script provides an example template construction
    and waveform cross-correlation clustering routine using most
    of the methods provided in this repository. 

    This particular example uses a set of "well-located" low
    frequency seismic events from Mt. Baker as the input 
    catalog events.

"""

# Provides a demonstration of an example workflow using `eqc_utils`
import os
from obspy import read_events
from obspy.clients.fdsn import Client
from eqcorrscan import Template, Tribe

from eqcorrscan_utils.logging import setup_terminal_logger
from eqcorrscan_utils.eventbank import initialize_event_bank
from eqcorrscan_utils.catalog import apply_phase_hints, filter_picks
from eqcorrscan_utils.template import augment_template, rename_templates
from eqcorrscan_utils.clustering import cluster
from eqcorrscan_utils.snuffler import plant

# Add 'snuffler' functionality to Tribe and Template class objects
plant()

# Set up command-line print-out for logging
Logger = setup_terminal_logger(__name__)

# Set up an IRIS webservices client
client = Client("IRIS")
Logger.info('')
# Load example event metadata from this repo
cat = read_events(os.path.join('.','mt_baker_catalog.xml'))

# Create an event bank & add catalog to it
ebank = initialize_event_bank(
    catalog=cat,
    base_path=os.path.join('..','demo_output','eventbank'),
    path_structure='{year}')                    

# Add phase hints from arrivals to picks
cat = apply_phase_hints(cat)

# Subset down to 3 nearest stations to Mt. Baker & only keep P picks
cat = filter_picks(
    catalog=cat,
    stations=['MBW','MBW2','SHUK'],
    phase_hints=['P'])

# Create a tribe of templates using IRIS webservices for waveforms
tribe = Tribe().construct(
    method='from_client',
    client_id=client,
    catalog=cat,
    samp_rate=50.,
    lowcut=0.5,
    highcut=20.,
    filt_order=4.,
    prepick=5.,
    length=40.,
    delayed=True,
    process_len=3600.,
    parallel=True,
    num_cores=6
)

# Rename the templates with EVID's including contributor codes
tribe = rename_templates(tribe)

# Augment the templates to have 3-C data (if applicable)
for template in tribe:
    template = augment_template(template, client)

# Visualize the tribe with snuffler
tribe.snuffle()

# Run clustering on templates
groups = cluster(
    tribe,
    corr_thresh=0.4,
    fill_value=1,
    show=True,
    save_subtribes=True,
    save_path=os.path.join('..','demo_output','template_clustering'))

