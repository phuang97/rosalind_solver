import logging

__author__ = "phuang97"
__email__ = "phuang0114@gmail.com"
__version__ = "0.0.1"

# without this logging will default to the logging.lastResort
# which has a level of WARNING and output stuff into sys.stderr for python3.2 +
logging.getLogger(__name__).addHandler(logging.NullHandler)
