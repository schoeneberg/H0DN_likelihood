"""Command-line interface and variant processing.

This module provides the CLI entry point for the h0_constrainer package, supporting
both single-configuration runs and batch variant analysis. It handles:
- Parsing command-line arguments
- Running the analysis with a single configuration file
- Running multiple variant analyses in isolated subprocesses
- Aggregating variant results into a single output file

Variant mode allows allows multiple analyses to be run by defining a base
configuration and multiple variant sections that override specific parameters.
Each variant runs in a fresh Python subprocess.

Command-line usage:
    h0_constrainer                    # Use config.ini in current directory
    h0_constrainer my_config.ini      # Use specified config file
    h0_constrainer -v variants.ini    # Run multiple variants

Functions:
    run_variant: Execute a single variant in a subprocess.
    cli_entry: Main CLI entry point for argument parsing and dispatch.
"""


import sys
import os
import tempfile
import configparser
import subprocess
from h0_constrainer.main import main  # used for the non-variant path
import h0_constrainer.config_reader as config_reader  

def run_variant(variant_name, base_config_dict, overrides, output_handle):
    """Execute a single variant analysis in an isolated subprocess.
    
    Creates a temporary configuration file by merging base config with variant
    overrides, then launches a fresh Python process to run the analysis. This
    is to ensure complete isolation between variants with no shared global state.
    
    Args:
        variant_name (str): Name of the variant (used for output labeling).
        base_config_dict (dict): Base configuration as nested dictionary
            (section -> {key: value}).
        overrides (dict): Parameter overrides for this variant (key: value pairs
            applied to DEFAULT section).
        output_handle (file): Open file handle for writing variant output.
    
    Side Effects:
        - Creates temporary .ini file (deleted after run)
        - Writes formatted output to output_handle including:
          * Variant name banner
          * Analysis stdout from subprocess
          * Stderr if non-zero exit code
          * Exit code if non-zero
    
    Note:
        Subprocess invokes: python -c "import h0_constrainer.cli as c; c.cli_entry()" <tmpfile>
        This ensures config_reader global state is fresh for each variant.
    """
    # Build a config from base + overrides
    config = configparser.ConfigParser()
    config.read_dict(base_config_dict)
    for key, value in overrides.items():
        config["DEFAULT"][key] = value

    # Write to a temporary .ini file
    with tempfile.NamedTemporaryFile("w", delete=False, suffix=".ini") as tmpfile:
        config_path = tmpfile.name
        config.write(tmpfile)

    try:
        # Write a banner for readability
        output_handle.write("=" * 80 + "\n")
        output_handle.write(f"Running variant: {variant_name}\n")
        output_handle.write("=" * 80 + "\n")

        # IMPORTANT:
        # Run a *fresh process* so no module/global state leaks between variants.
        #
        # We invoke the same CLI entrypoint by importing it in a -c command.
        # The config path is passed as an argv to that small command, which
        # cli_entry() will pick up via sys.argv[1:].
        proc = subprocess.run(
            [
                sys.executable,
                "-c",
                "import h0_constrainer.cli as c; c.cli_entry()",
                config_path,
            ],
            capture_output=True,
            text=True,
        )

        # Append the child process stdout/stderr to the variants_out.txt
        if proc.stdout:
            output_handle.write(proc.stdout)

        if proc.returncode != 0:
            # Make failures noisy and include stderr
            output_handle.write("\n[stderr]\n")
            output_handle.write(proc.stderr or "(no stderr captured)\n")
            output_handle.write(f"\n[exit code] {proc.returncode}\n")

        output_handle.write("\n")  # spacer after each variant

    finally:
        try:
            os.remove(config_path)
        except OSError:
            # If it couldn't be removed, don't crash the whole run
            pass


def cli_entry():
    """Main CLI entry point for h0_constrainer command.
    
    Parses command-line arguments and dispatches to either single-run mode
    or variant batch mode. This function is the entry point registered in
    pyproject.toml as the 'h0_constrainer' console script.
    
    Command-line modes:
        h0_constrainer                    # Single run with default config file name: config.ini
        h0_constrainer config.ini         # Single run with specified config
        h0_constrainer -v variants.ini    # Batch variant mode
    
    Variant File Format:
        [meta]
        base = config.ini
        
        [variant1]
        parameter1 = value1
        parameter2 = value2
        
        [variant2]
        parameter1 = value3
        ...
    
    Args:
        None (reads from sys.argv).
    
    Raises:
        ValueError: If variants.ini missing required [meta] section or base config.
        FileNotFoundError: If specified config files don't exist.
    
    Effects:
        - Single mode: Calls main() directly, output goes to stdout
        - Variant mode: Writes aggregated results to a named output file
                       derived from the variant file name.
    
    Example:
        $ h0_constrainer my_analysis.ini
        H₀=73.2±1.1 km/s/Mpc | χ²=145.2 | ndof=132 | ...
    """
    args = sys.argv[1:]

    # Variant mode: h0_constrainer -v variants.ini
    if len(args) >= 2 and args[0] == "-v":
        variant_file = args[1]

        vparser = configparser.ConfigParser()
        read_files = vparser.read(variant_file)
        if not read_files:
            raise FileNotFoundError(f"Variant file '{variant_file}' not found or unreadable.")

        if "meta" not in vparser or "base" not in vparser["meta"]:
            raise ValueError("variants.ini must include [meta] with base = config.ini")

        base_file = vparser["meta"]["base"]
        if not os.path.exists(base_file):
            raise FileNotFoundError(f"Base config file '{base_file}' not found.")

        base_parser = configparser.ConfigParser()
        read_bases = base_parser.read(base_file)
        if not read_bases:
            raise FileNotFoundError(f"Base config file '{base_file}' not found or unreadable.")

        # Build a dict snapshot of the base config (sections + DEFAULT)
        base_dict = {s: dict(base_parser.items(s)) for s in base_parser.sections()}
        base_dict["DEFAULT"] = dict(base_parser["DEFAULT"])


        # Derive output file name: basename without .ini + _out.txt
        variant_base = os.path.basename(variant_file)
        if variant_base.lower().endswith(".ini"):
            outname = variant_base[:-4] + "_out.txt"
        else:
            outname = variant_base + "_out.txt"

        # Run each variant in isolation; capture output to a single file
        with open(outname, "w", encoding="utf-8") as output:
            for section in vparser.sections():
                if section == "meta":
                    continue
                overrides = dict(vparser[section])
                run_variant(section, base_dict, overrides, output)

    else:
        # Non-variant path: run the tool once in this process
        config_file = args[0] if args else None
        main(config_file)