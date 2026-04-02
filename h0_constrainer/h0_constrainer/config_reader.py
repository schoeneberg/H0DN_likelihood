"""Configuration file reading and validation module.

This module provides utilities for loading and parsing INI-style configuration files,
tracking which configuration keys are read or skipped, and managing verbose/warning
output flags. It maintains global state for configuration tracking and output control.

Global Variables:
    read_keys (list): List of (key, value) tuples for configuration keys that were read.
    skipped_keys (list): List of (key, default) tuples for keys that were missing.
    verbose (bool): Flag controlling verbose output via vprint().
    warn (bool): Flag controlling warning output via wprint().
"""

import configparser

# Global lists to track which keys are read or skipped
read_keys = []
skipped_keys = []

def reset_keys():
    """Reset the global read/skipped key tracking lists.
    
    Call this at the start of each analysis run to ensure clean tracking
    when running multiple analyses in the same Python session.
    """
    global read_keys, skipped_keys
    read_keys = []
    skipped_keys = []

def load_config(config_file='config.ini'):
    """Load configuration from an INI file.
    
    Args:
        config_file (str, optional): Path to the configuration file. 
            Defaults to 'config.ini'.
    
    Returns:
        configparser.ConfigParser: Parsed configuration object.
    
    Example:
        >>> config = load_config('my_config.ini')
        >>> value = config['DEFAULT']['datadir']
    """
    config = configparser.ConfigParser()
    config.read(config_file)
    return config

def get_config_value(config, section, key, default=None, as_float=False, as_int=False, as_bool=False):
    """Retrieve a configuration value with automatic type conversion and tracking.
    
    This function attempts to read a configuration key from the specified section.
    If the key exists, it's added to the global read_keys list. If missing, the
    default value is used and added to the skipped_keys list.
    
    Args:
        config (configparser.ConfigParser): Configuration object.
        section (str): Section name in the INI file.
        key (str): Configuration key name.
        default (Any, optional): Default value if key is missing. Defaults to None.
        as_float (bool, optional): Convert value to float. Defaults to False.
        as_int (bool, optional): Convert value to int. Defaults to False.
        as_bool (bool, optional): Convert value to bool. Defaults to False.
    
    Returns:
        Any: Configuration value with appropriate type conversion, or default.
    
    Note:
        Only one of as_float, as_int, or as_bool should be True. If multiple are
        True, the priority is: as_float > as_int > as_bool.
    
    Example:
        >>> config = load_config()
        >>> vpec_error = get_config_value(config, 'DEFAULT', 'vpec_error', 240.0, as_float=True)
    """
    global read_keys, skipped_keys
    if key in config[section]:
        if as_float:
            value = config[section].getfloat(key)
        elif as_int:
            value = config[section].getint(key)
        elif as_bool:
            value = config[section].getboolean(key)
        else:
            value = config[section].get(key)
        read_keys.append((key, value))
        return value
    else:
        skipped_keys.append((key, default))
        return default

def print_config_keys():
    """Print summary of configuration keys that were read or skipped.
    
    Displays two lists:
    1. Keys that were successfully read from the configuration file.
    2. Keys that were missing and used default values.
    
    This is useful for debugging configuration issues and understanding
    which parameters are active in a run.
    """
    print("\n **Config keys that were read:**")
    for key, value in read_keys:
        print(f"  - {key}: {value}")
    print("\n **Config keys that were missing or skipped (using default values):**")
    for key, default in skipped_keys:
        print(f"  - {key}: (Default: {default})")

def give_keys():
    """Return the lists of read and skipped configuration keys.
    
    Returns:
        tuple: (read_keys, skipped_keys) where each is a list of tuples.
            - read_keys: List of (key, value) tuples.
            - skipped_keys: List of (key, default) tuples.
    """
    return read_keys, skipped_keys


def set_flags(*, verbose_flag=None, warn_flag=None):
    """Set global verbose and warning output flags.
    
    This function should be called once from main() after reading the configuration
    to initialize the global output control flags.
    
    Args:
        verbose_flag (bool, optional): Enable verbose output. If None, flag unchanged.
        warn_flag (bool, optional): Enable warning output. If None, flag unchanged.
    
    Example:
        >>> set_flags(verbose_flag=True, warn_flag=True)
        >>> vprint("This will be printed")
    """
    global verbose, warn
    if verbose_flag is not None:
        verbose = bool(verbose_flag)
    if warn_flag is not None:
        warn = bool(warn_flag)


def vprint(*args, **kwargs):
    """Print only when verbose mode is enabled.
    
    Wrapper around print() that only outputs when the global verbose flag is True.
    Accepts all arguments that print() accepts.
    
    Args:
        *args: Positional arguments passed to print().
        **kwargs: Keyword arguments passed to print().
    """
    if verbose:
        print(*args, **kwargs)

def wprint(*args, **kwargs):
    """Print only when warning mode is enabled.
    
    Wrapper around print() that only outputs when the global warn flag is True.
    Accepts all arguments that print() accepts. Used for warnings that may not
    be critical but should be visible to users who enable warnings.
    
    Args:
        *args: Positional arguments passed to print().
        **kwargs: Keyword arguments passed to print().
    """
    if warn:
        print(*args, **kwargs)


