import logging

import pandas as pd

from obspy import Inventory
from pyrocko import cake

from eqcutil.util.logging import rich_error_message
from eqcutil.catalog.model_phases import make_model, model_raypaths, ray_summary

class RayModeler(object):
    """
    A class for managing the inputs and outputs of modeling 
    """
    def __init__(self, model='P4', inventory=None):
        # Set up logging
        self.Logger = logging.getLogger(self.__name__())
        # Compatability checks on model
        if isinstance(model, cake.LayeredModel):
            self.model = model
        elif isinstance(model, str):
            try:
                self.make_model(name=model)
            except Exception as e:
                self.Logger.critical(rich_error_message(e))
        else:
            self.Logger.critical(f'TypeError: "model" value {model} not supported.')
        
        # inventory compatability checks
        if isinstance(inventory, Inventory):
            self.inv = inventory
            if len(self.inv) == 0:
                self.Logger.warning('Initialized with an empty inventory - will need to add stations to model rays!')
        elif inventory is None:
            self.Logger.warning('Initialized without an inventory - will need to add an inventory to model rays!')
        else:
            self.Logger.critical(f'TypeError: "inventory" type {type(inventory)} not supported.')
        # Initialize Summary Object
        self.summary = pd.DataFrame()
        # Initialize origin list
        self.origins = []

    def __name__(self):
        return 'eqcorrscan_utils.cake.RayModeler'

    def make_model(self, name='P4'):
        try:
            self.model = make_model(name=name)
        except Exception as e:
            raise

    def model_raypaths(self, origin, phases=['P','S']):
        if not (origin, phases) in self.origins:
            results = model_raypaths(self.model, origin, self.inv, phases=phases)
            results = ray_summary(results, origin, self.inv)
            self.summary = pd.concat([self.summary, results], ignore_index=True)
            self.origins.append((origin, phases))
        else:
            self.Logger.warning('This origin-phases combination has already been run!')
    
