# Provides a demonstration of an example workflow using `eqc_utils`

from obspy.clients.fdsn import Client
from eqcorrscan import Tribe

from eqc_utils.

# Connect to a USGS webservices client (lets us get pick info from COMCAT)
client = Client("USGS")

# Get a few catalog events from the FDSN

