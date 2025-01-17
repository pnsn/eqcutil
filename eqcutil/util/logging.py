"""
:module: eqcorrscan_utils.logging
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

import sys, logging

class CriticalExitHandler(logging.Handler):
    """A custom :class:`~logging.Handler` sub-class that emits a sys.exit
    if a logging instance emits a logging.CRITICAL level message

    Constructed through a prompt to ChatGPT and independently
    tested by the author.

    """
    def __init__(self, exit_code=1, **kwargs):
        super().__init__(**kwargs)
        self.exit_code = exit_code
        
    def emit(self, record):
        if record.levelno == logging.CRITICAL:
            sys.exit(self.exit_code)

def rich_error_message(e, additional_text=None):
    """Given the raw output of an "except Exception as e"
    formatted clause, return a string that emulates
    the error message print-out from python

    e.g., 
    'TypeError: input "input" is not type int'

    :param e: Error object
    :type e: _type_
    :return: error message
    :rtype: str
    """    
    etype = type(e).__name__
    emsg = str(e)
    if additional_text is not None:
        return f'{etype}: {emsg} {additional_text}'
    else:
        return f'{etype}: {emsg}'


def basic_logger_config(level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"):
    logging.basicConfig(level=level, format=format)


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

