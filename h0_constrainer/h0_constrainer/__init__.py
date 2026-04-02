"""H0DN: Distance network analysis package for constraining the Hubble constant (H0).

This package is a Python implementation developed for the H0DN collaboration
based on the original IDL code written by Stefano Casertano.

Modules:
    config_reader: Configuration file parsing and validation
    data_loader: Data file loading and preprocessing
    equations: Linear equation system construction
    solver: Weighted least-squares solver
    intercept: Hubble flow alpha intercept computation
    logger: Results logging and serialization
    main: Primary analysis workflow
    cli: Command-line interface

Usage:
    Command-line:
        $ h0_constrainer config.ini
        $ h0_constrainer -v variants.ini
    
    Programmatic:
        >>> from h0_constrainer import main
        >>> main.main('my_config.ini')

Version: 1.0
Author: [Project team]
"""

from . import config_reader
from . import data_loader
from . import equations
from . import intercept
from . import solver
from . import logger
from . import main
from . import cli