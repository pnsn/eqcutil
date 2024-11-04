"""
:module: eqc_utils.logging
:auth: Nathan T. Stevens
:email: ntsteven@uw.edu
:org: Pacific Northwest Seismic Network
:license: GNU GPLv3
:purpose: 
    This module contains helper methods for setting up logging
    to the command line

:attribution: 
    These method(s) are abstracted from the Python Logging Cookbook
    https://docs.python.org/3/howto/logging-cookbook.html
    And use a bug-fix provided by StackExchange user "Euclides"
    from Sep 8, 2015 for preventing duplicate StreamHandler
    instances.
"""

import os, logging

def setup_terminal_logger(name, level=logging.INFO):
    """QuickStart setup for a write-to-command-line-only
    logger that has safety catches against creating
    multiple StreamHandlers with repeat calls of a
    script in ipython or similar interactive python session.

    :param name: name for the logger
        e.g., __name__ is fairly customary
    :type name: str
    :param level: logging level, defaults to logging.INFO
    :type level: int, optional
    :return: Logger
    :rtype: logging.RootLogger
    """    
    ## SET UP LOGGING
    Logger = logging.getLogger(name)
    # Set logging level to INFO
    Logger.setLevel(level)

    # Prevent duplication during testing
    # Solution from https://stackoverflow.com/questions/31403679/python-logging-module-duplicated-console-output-ipython-notebook-qtconsole
    # User Euclides (Sep 8, 2015)
    handler_console = None
    handlers = Logger.handlers
    for h in handlers:
        if isinstance(h, logging.StreamHandler):
            handler_console = h
            break
    # Set up logging to terminal
    if handler_console is None:
        ch = logging.StreamHandler()
        # Set up logging line format
        fmt = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        ch.setFormatter(fmt)
        # Add formatting & handler
        Logger.addHandler(ch)
    return Logger

