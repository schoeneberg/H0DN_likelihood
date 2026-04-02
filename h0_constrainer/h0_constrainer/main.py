"""Main entry point and workflow.

This module coordinates the complete distance network analysis workflow:
1. Load and parse configuration file
2. Load and filter data files (hosts, calibrators, Hubble flow, etc.)
3. Apply exclusion/inclusion criteria provided in the configuration file
4. Build the equation system with all constraints
5. Solve for H0, absolute magnitudes, intercepts, and other parameters
6. Optionally save results

The analysis framework is designed to incorporate a comprehensive distance network:
- Geometric distance measurements to anchors
- Host galaxy distances from various primary distance indicators (Cepheids, TRGB, JAGB, Miras)
- Calibrator magnitude data for secondary indicators (SNe Ia, SBF, TF, SNe II)
- Hubble flow measurements for inferring H0 from absolute magnitudes
- Additional constraints, such as megamasers, EPM, Coma cluster, and galaxy group data
- Proper tracking of shared uncertainties/covariances between different measurements

Command-line usage:
    h0_constrainer [config_file]
    h0_constrainer -v variants.ini

Functions:
    main: Primary analysis workflow.
"""

import os
import sys
import numpy as np
import pandas as pd
from h0_constrainer import config_reader, data_loader, equations, solver, logger




def main(cf_file = None):
    """Execute the complete distance ladder analysis workflow.
    
    Orchestrates all steps from configuration loading through result output:
    1. Configuration: Load and validate configuration parameters
    2. Data Loading: Load all input data files (hosts, calibrators, Hubble flow, etc.)
    3. Data Filtering: Apply inclusion/exclusion criteria and filter different types of data to ensure consistency
    4. Equation Building: Construct linear system with all constraints and covariances
    5. Solution: Solve generalized least-squares system for H0 and parameters
    6. Optional output: Save results to log file, JSON
    
    Args:
        cf_file (str, optional): Path to configuration file. If None (default), uses command-line
            argument or defaults to 'config.ini'. 
    
    Configuration File Requirements:
        Extensive examples of the configuration file are provided in ../h0_constrainer/configs/config.ini

        Parameters:
        - datadir: Directory containing data files
        - host_data: Host galaxy distance measurements file
        - Calibrators (sn1a_calib, sn2_calib, tf_calib, or sbf_calib)
        - Corresponding Hubble flow file or alpha parameters for secondary indicators
        - verbose, warn: Output verbosity flags
        - anchors: Geometric anchor measurements file
        - groups: Galaxy group membership file
        - mm, epm, coma: Additional geometric constraints
        - host_exclude_list, calib_exclude: Exclusion patterns
        - host_include_list, calib_include: Inclusion patterns
        - q0, j0: Cosmographic parameters
        - vpec_error: Peculiar velocity uncertainty
        - zmin_hf, zmax_hf: Redshift range for Hubble flow
        - logfile, savefile: Output file paths
        - Many method-specific parameters (see configuration examples in ../h0_constrainer/configs/config.ini)
    
    Raises:
        ValueError: If configuration is invalid (e.g., calibrator without HF/alpha,
            HF/alpha without calibrator, missing required files).
        FileNotFoundError: If specified data files don't exist.
    
    Effects:
        - Prints results summary to stdout (format depends on verbose setting)
        - Writes log file if logfile specified
        - Saves detailed results JSON if savefile specified
        - Updates global config_reader state (read_keys, skipped_keys, verbose, warn)
    
    Example:
        >>> main('my_analysis.ini')
        H₀=73.2±1.1 km/s/Mpc | χ²=145.2 | ndof=132 | χ²_red=1.10 | M_B=-19.25±0.03
    
    Note:
        When called via command-line interface (h0_constrainer), sys.argv is used
        to determine the config file. When called programmatically, cf_file can be
        specified directly.
    """
    if cf_file:
      config_file = cf_file
    elif len(sys.argv)==2:
      config_file = sys.argv[1]
      if not os.path.exists(config_file):
        raise ValueError(f"You passed the config file '{config_file}' as a command line argument, but I could not find the specified file!")
    else:
      config_file = "config.ini"

    # Reset key tracking for clean state (important for batch/variant runs)
    config_reader.reset_keys()

    # Load configuration
    config = config_reader.load_config(config_file)
    datadir = config_reader.get_config_value(config, "DEFAULT", "datadir", "./")
    verbose = config_reader.get_config_value(config, "DEFAULT", "verbose", False, as_bool=True)
    warn = config_reader.get_config_value(config, "DEFAULT", "warn", True, as_bool=True)

    config_reader.set_flags(verbose_flag=verbose, warn_flag=warn)

    anchor_file = config_reader.get_config_value(config, "DEFAULT", "anchors", "anchors.dat").strip().lower()
    host_file = config_reader.get_config_value(config, "DEFAULT", "host_data", "").strip().lower()

    sn1a_calib_file = config_reader.get_config_value(config, "DEFAULT", "sn1a_calib", "").strip().lower()
    sn2_calib_file = config_reader.get_config_value(config, "DEFAULT", "sn2_calib", "").strip().lower()
    tf_calib_file = config_reader.get_config_value(config, "DEFAULT", "tf_calib", "").strip().lower()
    sbf_calib_file = config_reader.get_config_value(config, "DEFAULT", "sbf_calib", "").strip().lower()

    groups_file = config_reader.get_config_value(config, "DEFAULT", "groups", "").strip().lower()

    sn1a_hf_file = config_reader.get_config_value(config, "DEFAULT", "sn1a_hf", "").strip().lower()
    sn1a_hf_cov_file = config_reader.get_config_value(config, "DEFAULT", "sn1a_hf_cov", "").strip().lower()
    sn2_hf_file = config_reader.get_config_value(config, "DEFAULT", "sn2_hf", "").strip().lower()
    sn2_hf_cov_file = config_reader.get_config_value(config, "DEFAULT", "sn2_hf_cov", "").strip().lower()
    sn1a_ignore_offdiag = config_reader.get_config_value(config, "DEFAULT", "sn1a_ignore_offdiag", False, as_bool=True)
    sn1a_IR = config_reader.get_config_value(config, "DEFAULT", "sn1a_IR", False, as_bool=True)
    sn1a_intrinsic = config_reader.get_config_value(config, "DEFAULT", "sn1a_intrinsic", 0.0, as_float=True)
    tf_hf_file = config_reader.get_config_value(config, "DEFAULT", "tf_hf", "").strip().lower()
    tf_hf_cov_file = config_reader.get_config_value(config, "DEFAULT", "tf_hf_cov", "").strip().lower()
    sbf_hf_file = config_reader.get_config_value(config, "DEFAULT", "sbf_hf", "").strip().lower()
    sbf_hf_cov_file = config_reader.get_config_value(config, "DEFAULT", "sbf_hf_cov", "").strip().lower()

    mm_file = config_reader.get_config_value(config, "DEFAULT", "mm", "").strip().lower()
    epm_file = config_reader.get_config_value(config, "DEFAULT", "epm", "").strip().lower()    
    coma_file = config_reader.get_config_value(config, "DEFAULT", "coma", "").strip().lower()
    sigma_mag_fp = config_reader.get_config_value(config, "DEFAULT", "sigma_mag_fp", 0.0371, as_float=True)

    mm_dz = config_reader.get_config_value(config, "DEFAULT", "mm_dz", True, as_bool=True)
    epm_dz = config_reader.get_config_value(config, "DEFAULT", "epm_dz", True, as_bool=True)

    logfile = config_reader.get_config_value(config, "DEFAULT", "logfile", "none").strip()
    savefile = config_reader.get_config_value(config, "DEFAULT", "savefile", "none").strip()
    # Exclusion lists
    host_exclude_list_str = config_reader.get_config_value(config, "DEFAULT", "host_exclude_list", "")
    host_exclude_list = [item.strip() for item in host_exclude_list_str.split(",")] if host_exclude_list_str else []
    calib_exclude_str = config_reader.get_config_value(config, "DEFAULT", "calib_exclude", "")
    calib_exclude = [item.strip() for item in calib_exclude_str.split(",")] if calib_exclude_str else []
    # Inclusion lists
    host_include_list_str = config_reader.get_config_value(config, "DEFAULT", "host_include_list", "")
    host_include_list = [item.strip() for item in host_include_list_str.split(",")] if host_include_list_str else []
    calib_include_str = config_reader.get_config_value(config, "DEFAULT", "calib_include", "")
    calib_include = [item.strip() for item in calib_include_str.split(",")] if calib_include_str else []

    q0 = config_reader.get_config_value(config, "DEFAULT", "q0", -0.55, as_float=True)
    j0 = config_reader.get_config_value(config, "DEFAULT", "j0", 1.0, as_float=True)

    zmin_hf = config_reader.get_config_value(config, "DEFAULT", "zmin_hf", 0.0, as_float=True)
    zmax_hf = config_reader.get_config_value(config, "DEFAULT", "zmax_hf", 100.0, as_float=True)
    if zmin_hf == 0.0 and zmax_hf == 100.0:
        zbounds_hf = None
    else:
        zbounds_hf = (zmin_hf, zmax_hf)
    z_range_sn1a = config_reader.get_config_value(config, "DEFAULT", "z_range_sn1a", False, as_bool=True)
    z_range_sn2 = config_reader.get_config_value(config, "DEFAULT", "z_range_sn2", False, as_bool=True)
    z_range_tf = config_reader.get_config_value(config, "DEFAULT", "z_range_tf", False, as_bool=True)
    z_range_sbf = config_reader.get_config_value(config, "DEFAULT", "z_range_sbf", False, as_bool=True)

    alpha_sn1a = config_reader.get_config_value(config, "DEFAULT", "alpha_sn1a", "", as_float=True)
    alpha_sn1a_error = config_reader.get_config_value(config, "DEFAULT", "alpha_sn1a_error", "", as_float=True)
    alpha_sn2 = config_reader.get_config_value(config, "DEFAULT", "alpha_sn2", "", as_float=True)
    alpha_sn2_error = config_reader.get_config_value(config, "DEFAULT", "alpha_sn2_error", "", as_float=True)
    alpha_tf = config_reader.get_config_value(config, "DEFAULT", "alpha_tf", "", as_float=True)
    alpha_tf_error = config_reader.get_config_value(config, "DEFAULT", "alpha_tf_error", "", as_float=True)
    alpha_sbf = config_reader.get_config_value(config, "DEFAULT", "alpha_sbf", "", as_float=True)
    alpha_sbf_error = config_reader.get_config_value(config, "DEFAULT", "alpha_sbf_error", "", as_float=True)

    unused_hosts_via_covariance = config_reader.get_config_value(config, "DEFAULT", "unused_hosts_via_covariance", False, as_bool=True)

    vcorr_mm = config_reader.get_config_value(config, "DEFAULT", "vcorr_mm", "2M++")
    vcorr_epm = config_reader.get_config_value(config, "DEFAULT", "vcorr_epm", "2M++")
    vpec_error = config_reader.get_config_value(config, "DEFAULT", "vpec_error", 240.0, as_float=True)
    sbf_vpec_error = config_reader.get_config_value(config, "DEFAULT", "sbf_vpec_error", vpec_error, as_float=True)
    sn1a_vp = config_reader.get_config_value(config, "DEFAULT", "sn1a_vp", "2m++")
    tf_v = config_reader.get_config_value(config, "DEFAULT", "tf_v", "2m++")
    sbf_v = config_reader.get_config_value(config, "DEFAULT", "sbf_v", "2m++")

    if verbose:
        config_reader.print_config_keys()

    read, skipped = config_reader.give_keys()
    read = [r[0] for r in read]
    #TODO THIS PART IS TO BE WRITTEN
    # if "rung3" not in read and any([x in read for x in ["a_b","q0","j0","sigma_a_b","rung_hf"]]):
    #   raise ValueError("If rung3 is not provided, you cannot provide any of 'a_b','sigma_a_b','q0','j0','rung_hf'")
    # if "rung_hf" in read and any([x in read for x in ["a_b","sigma_a_b"]]):
    #   raise ValueError("If rung_hf is provided, you cannot provide any of 'a_b','sigma_a_b'")
    # if ("a_b" in read or "sigma_a_b" in read) and any([x in read for x in ["rung_hf","q0","j0"]]):
    #   raise ValueError("If a_b or sigma_a_b is provided, you cannot provide any of 'rung_hf','q0','j0'")
  
    # Validate that no HF or alpha values are provided without their corresponding calibration file
    if "sn1a_calib" not in read and any(k in read for k in ["sn1a_hf", "alpha_sn1a", "alpha_sn1a_error"]): 
        raise ValueError("If 'sn1a_calib' file is not provided, you cannot provide 'sn1a_hf' or ('alpha_sn1a' and 'alpha_sn1a_error')")
    if "sn2_calib" not in read and any(k in read for k in ["sn2_hf", "alpha_sn2", "alpha_sn2_error"]): 
        raise ValueError("If 'sn2_calib' file is not provided, you cannot provide 'sn2_hf' or ('alpha_sn2' and 'alpha_sn2_error')")
    if "tf_calib" not in read and any(k in read for k in ["tf_hf", "alpha_tf", "alpha_tf_error"]): 
        raise ValueError("If 'tf_calib' file is not provided, you cannot provide 'tf_hf' or ('alpha_tf' and 'alpha_tf_error')")
    if "sbf_calib" not in read and any(k in read for k in ["sbf_hf", "alpha_sbf", "alpha_sbf_error"]):
        raise ValueError("If 'sbf_calib' file is not provided, you cannot provide 'sbf_hf' or ('alpha_sbf' and 'alpha_sbf_error')")
    
    # Validate that each calibration file is supported by either a hubble flow file or alpha+error
    # TODO: For Coma cluster analysis, it may be necessary to provide sn1a_calib without HF or alpha info
    if "sn1a_calib" in read and not ("sn1a_hf" in read or ("alpha_sn1a" in read and "alpha_sn1a_error" in read)):
        raise ValueError("If 'sn1a_calib' file is provided, you must also provide either 'sn1a_hf' or both 'alpha_sn1a' and 'alpha_sn1a_error'")
    if "sn2_calib" in read and not ("sn2_hf" in read or ("alpha_sn2" in read and "alpha_sn2_error" in read)): 
        raise ValueError("If 'sn2_calib' file is provided, you must also provide either 'sn2_hf' or both 'alpha_sn2' and 'alpha_sn2_error'")
    if "tf_calib" in read and not ("tf_hf" in read or ("alpha_tf" in read and "alpha_tf_error" in read)): 
        raise ValueError("If 'tf_calib' file is provided, you must also provide either 'tf_hf' or both 'alpha_tf' and 'alpha_tf_error'")
    if "sbf_calib" in read and not ("sbf_hf" in read or ("alpha_sbf" in read and "alpha_sbf_error" in read)):
        raise ValueError("If 'sbf_calib' file is provided, you must also provide either 'sbf_hf' or both 'alpha_sbf' and 'alpha_sbf_error'")



    # Initialize data containers
    host_df = None
    sn1a_calib_df = None
    sn2_calib_df = None
    tf_calib_df = None
    sbf_calib_df = None
    sn1a_hf_df = None
    sn1a_hf_cov = None
    sn2_hf_df = None
    sn2_hf_cov = None
    tf_hf_df = None
    tf_hf_cov = None
    sbf_hf_df = None
    sbf_hf_cov = None

    # Load data files
    host_df = data_loader.load_hosts(datadir, host_file)
    if sn1a_IR:
        sn1a_calib_df = data_loader.load_calibrators(datadir, sn1a_calib_file, "sn_IR")
    else:
        sn1a_calib_df = data_loader.load_calibrators(datadir, sn1a_calib_file, "sn")
    if sn1a_calib_df is not None:
        sn1a_calib_df["sigma"] = np.sqrt(sn1a_calib_df["sigma"]**2 + sn1a_intrinsic**2)

    sn2_calib_df = data_loader.load_calibrators(datadir, sn2_calib_file, "sn")
    tf_calib_df = data_loader.load_calibrators(datadir, tf_calib_file, "tf")
    sbf_calib_df = data_loader.load_calibrators(datadir, sbf_calib_file, "sbf")

    groups_df = data_loader.load_groups(datadir, groups_file)

    # Host data filtering (apply inclusion/exclusion lists)
    if host_df is not None:
        if host_exclude_list and host_include_list:
            config_reader.wprint("\n=== WARNING === Both inclusion and exclusion lists are provided for host_data. Ignoring inclusion list.")
            host_df = host_df[~host_df.apply(
                lambda row: data_loader.should_exclude(row, host_exclude_list, ["host", "method", "anchor", "source"]),
                axis=1)]
        elif host_exclude_list:
            host_df = host_df[~host_df.apply(
                lambda row: data_loader.should_exclude(row, host_exclude_list, ["host", "method", "anchor", "source"]),
                axis=1)]
        elif host_include_list:
            host_df = host_df[host_df.apply(
                lambda row: data_loader.should_include(row, host_include_list, ["host", "method", "anchor", "source"]),
                axis=1)]

    # SNe Ia calibrator filtering
    if sn1a_calib_df is not None:
        if calib_exclude and calib_include:
            config_reader.wprint("\n=== WARNING === Both inclusion and exclusion lists are provided for sn1a_calib. Ignoring inclusion list.")
            sn1a_calib_df = sn1a_calib_df[~sn1a_calib_df.apply(
                lambda row: data_loader.should_exclude(row, calib_exclude, ["host", "SN"]),
                axis=1)]
        elif calib_exclude:
            sn1a_calib_df = sn1a_calib_df[~sn1a_calib_df.apply(
                lambda row: data_loader.should_exclude(row, calib_exclude, ["host", "SN"]),
                axis=1)]
        elif calib_include:
            sn1a_calib_df = sn1a_calib_df[sn1a_calib_df.apply(
                lambda row: data_loader.should_include(row, calib_include, ["host", "SN"]),
                axis=1)]
    
    # SNe II calibrator filtering
    if sn2_calib_df is not None:
        if calib_exclude and calib_include:
            config_reader.wprint("\n=== WARNING === Both inclusion and exclusion lists are provided for sn2_calib. Ignoring inclusion list.")
            sn2_calib_df = sn2_calib_df[~sn2_calib_df.apply(
                lambda row: data_loader.should_exclude(row, calib_exclude, ["host", "SN"]),
                axis=1)]
        elif calib_exclude:
            sn2_calib_df = sn2_calib_df[~sn2_calib_df.apply(
                lambda row: data_loader.should_exclude(row, calib_exclude, ["host", "SN"]),
                axis=1)]
        elif calib_include:
            sn2_calib_df = sn2_calib_df[sn2_calib_df.apply(
                lambda row: data_loader.should_include(row, calib_include, ["host", "SN"]),
                axis=1)]

    # Tully-Fisher calibrator filtering
    if tf_calib_df is not None:
        if calib_exclude and calib_include:
            config_reader.wprint("\n=== WARNING === Both inclusion and exclusion lists are provided for tf_calib. Ignoring inclusion list.")
            tf_calib_df = tf_calib_df[~tf_calib_df.apply(
                lambda row: data_loader.should_exclude(row, calib_exclude, ["host", "PGC", "anchor", "source", "method"]),
                axis=1)]
        elif calib_exclude:
            tf_calib_df = tf_calib_df[~tf_calib_df.apply(
                lambda row: data_loader.should_exclude(row, calib_exclude, ["host", "PGC", "anchor", "source", "method"]),
                axis=1)]
        elif calib_include:
            tf_calib_df = tf_calib_df[tf_calib_df.apply(
                lambda row: data_loader.should_include(row, calib_include, ["host", "PGC", "anchor", "source", "method"]),
                axis=1)]

    # SBF calibrator filtering
    if sbf_calib_df is not None:
        if calib_exclude and calib_include:
            config_reader.wprint("\n=== WARNING === Both inclusion and exclusion lists are provided for sbf_calib. Ignoring inclusion list.")
            sbf_calib_df = sbf_calib_df[~sbf_calib_df.apply(
                lambda row: data_loader.should_exclude(row, calib_exclude, ["host"]),
                axis=1)]
        elif calib_exclude:
            sbf_calib_df = sbf_calib_df[~sbf_calib_df.apply(
                lambda row: data_loader.should_exclude(row, calib_exclude, ["host"]),
                axis=1)]
        elif calib_include:
            sbf_calib_df = sbf_calib_df[sbf_calib_df.apply(
                lambda row: data_loader.should_include(row, calib_include, ["host"]),
                axis=1)]

    # Display filtered data counts
    if verbose:
        if host_df is not None:
            print(f"\nHosts:  {len(host_df)} rows after filtering exclusions/inclusions.")
        if groups_df is not None:
            print(f"Groups:  {len(groups_df)} rows after filtering exclusions/inclusions.")
        if sn1a_calib_df is not None:
            print(f"SNe Ia: {len(sn1a_calib_df)} rows after filtering exclusions/inclusions.")
        if sn2_calib_df is not None:
            print(f"SNe II: {len(sn2_calib_df)} rows after filtering exclusions/inclusions.")
        if tf_calib_df is not None:
            print(f"TF:     {len(tf_calib_df)} rows after filtering exclusions/inclusions.")
        if sbf_calib_df is not None:
            print(f"SBF:    {len(sbf_calib_df)} rows after filtering exclusions/inclusions.")

    # Filter groups to only include those with hosts in the data
    if host_df is not None:
        if groups_df is not None:
            groups_df = data_loader.filter_groups(host_df, groups_df, sn1a_calib_df=sn1a_calib_df, sn2_calib_df=sn2_calib_df, tf_calib_df=tf_calib_df, sbf_calib_df=sbf_calib_df)
    
    # Initialize removed entries lists
    removed_hosts_without_calibrators=None,
    removed_sn1a_calib_without_hosts=None,
    removed_sn2_calib_without_hosts=None,
    removed_tf_calib_without_hosts=None,
    removed_sbf_calib_without_hosts=None,

    config_reader.vprint("\n", end="")
    # Drop calibrators for hosts that are absent in host data
    if sn1a_calib_df is not None:
        config_reader.vprint("SNe Ia", end=" ")
        sn1a_calib_df, removed_sn1a_calib_without_hosts = data_loader.rung_equalize(host_df, sn1a_calib_df, df_groups=groups_df)
    if sn2_calib_df is not None:
        config_reader.vprint("SNe II", end=" ")
        sn2_calib_df, removed_sn2_calib_without_hosts = data_loader.rung_equalize(host_df, sn2_calib_df, df_groups=groups_df)
    if tf_calib_df is not None:
        config_reader.vprint("TF    ", end=" ")
        tf_calib_df, removed_tf_calib_without_hosts = data_loader.rung_equalize(host_df, tf_calib_df, df_groups=groups_df)
    if sbf_calib_df is not None:
        config_reader.vprint("SBF   ", end=" ")
        sbf_calib_df, removed_sbf_calib_without_hosts = data_loader.rung_equalize(host_df, sbf_calib_df, df_groups=groups_df)

    # Drop hosts that have no corresponding calibrators (unless unused_hosts_via_covariance is True)
    if not unused_hosts_via_covariance:
        combined_hosts = set()
        if sn1a_calib_df is not None:
            combined_hosts.update(sn1a_calib_df["host"])
        if sn2_calib_df is not None:
            combined_hosts.update(sn2_calib_df["host"])
        if tf_calib_df is not None:
            combined_hosts.update(tf_calib_df["host"])
        if sbf_calib_df is not None:
            combined_hosts.update(sbf_calib_df["host"])
        combined_df = pd.DataFrame({"host": list(combined_hosts)})
        host_df, removed_hosts_without_calibrators = data_loader.rung_equalize(host_df, combined_df, reverse=True, df_groups=groups_df)

    # Drop anchor rows and obtain MAS (Method-Anchor-Source) reference arrays
    (host_df,
    mas_ref_error, mas_ref_value, mas_ref_anchor,
    mas_ref_method, mas_ref_source, n_mas, mas, mas_index, mas_error) = data_loader.drop_anchor_rows(host_df)

    
    # Load anchors data
    anchor_data, anchor_names, anchor_g_mu_sigma = data_loader.load_anchors(datadir, anchorfile=anchor_file)


    # --- SN Ia load Hubble flow ---
    sn1a_vp_column = None
    if sn1a_hf_file:
        if sn1a_hf_cov_file:
            sn1a_hf_df, sn1a_hf_cov = data_loader.load_hf(
                datadir, sn1a_hf_file, "sn1a", covar_filename=sn1a_hf_cov_file
            )
        else:
            method = "sn1a_IR" if sn1a_IR else "sn1a"
            sn1a_hf_df, sn1a_hf_cov = data_loader.load_hf(datadir, sn1a_hf_file, method)
        sn1a_hf_df["mb_err"] = np.sqrt(sn1a_hf_df["mb_err"]**2 + sn1a_intrinsic**2)

        # Usual SN Ia options (non-IR mode)
        allowed_sn1a_usual = ["plain", "2m++", "2mrs", "sdss", "cmb"]
        sn1a_map_usual = {
            "plain": "vp",
            "2m++":  "vp_2mpp",
            "2mrs":  "vp_2mrs",
            "sdss":  "vp_2mpp_sdss_6df",
            "cmb":   "vp_cmb",
        }

        # IR mode options (restricted)
        allowed_sn1a_ir = ["2m++", "cmb"]
        sn1a_map_ir = {
            "2m++": "vp_2mpp",
            "cmb":  "vp_cmb",
        }

        key = sn1a_vp.strip().lower()

        if not sn1a_IR:
            # ---- Non-IR: use full SN Ia options (including 'plain') ----
            if key not in allowed_sn1a_usual:
                raise ValueError(
                    f"Invalid choice for sn1a_vp '{sn1a_vp}'. Must be one of: {allowed_sn1a_usual}"
                )
            sn1a_vp_column = sn1a_map_usual[key]
            # Special rule: if 'cmb' is requested but the column is absent, create it as zeros
            if sn1a_vp_column == "vp_cmb" and sn1a_vp_column not in sn1a_hf_df.columns:
                sn1a_hf_df[sn1a_vp_column] = 0.0

            if sn1a_vp_column not in sn1a_hf_df.columns:
                raise KeyError(
                    f"Column '{sn1a_vp_column}' not found in SN Ia Hubble-flow file. "
                    f"Available columns include: {list(sn1a_hf_df.columns)}"
            )
        else:
            # ---- IR mode: only '2m++' and 'cmb' are allowed ----
            if key in allowed_sn1a_usual and key not in allowed_sn1a_ir:
                # user picked a normally valid option that's not allowed in IR
                raise ValueError("For SN Ia IR mode, only '2m++' and 'CMB' are possible.")
            if key not in allowed_sn1a_ir:
                raise ValueError(
                    f"Invalid choice for sn1a_vp '{sn1a_vp}'. For SN Ia IR mode, "
                    f"must be one of: {allowed_sn1a_ir}"
                )
            sn1a_vp_column = sn1a_map_ir[key]




    # --- SN II load Hubble flow ---
    if sn2_hf_file:
        sn2_hf_df, sn2_hf_cov = data_loader.load_hf(
            datadir, sn2_hf_file, "sn2", covar_filename=sn2_hf_cov_file
        )

    # --- Tully–Fisher load Hubble flow ---
    tf_v_column = None
    if tf_hf_file:
        tf_hf_df, tf_hf_cov = data_loader.load_hf(
            datadir, tf_hf_file, "tf", covar_filename=tf_hf_cov_file
        )

        allowed_tf = ["cmb", "2mrs", "2m++", "sdss"]
        tf_map = {
            "cmb":  "v_cmb",
            "2mrs": "v_2mrs",
            "2m++": "v_2mpp",
            "sdss": "v_2mpp_sdss_6df",
        }

        key = tf_v.strip().lower()
        if key not in allowed_tf:
            raise ValueError(
                f"Invalid choice for tf_v '{tf_v}'. Must be one of: {allowed_tf}"
            )

        tf_v_column = tf_map[key]
        if tf_v_column not in tf_hf_df.columns:
            raise KeyError(
                f"Column '{tf_v_column}' not found in TF Hubble-flow file. "
                f"Available columns include: {list(tf_hf_df.columns)}"
            )


    # --- SBF load Hubble flow ---
    sbf_v_column = None
    if sbf_hf_file:
        sbf_hf_df, sbf_hf_cov = data_loader.load_hf(
            datadir, sbf_hf_file, "sbf", covar_filename=sbf_hf_cov_file
        )

        # SBF supports only: 'cmb' (mapped to group frame) and '2m++'
        allowed_sbf = ["cmb", "2m++"]
        sbf_map = {
            "cmb":  "vgrp",   # 'cmb' input -> group-frame column
            "2m++": "v2m++",
        }

        key = sbf_v.strip().lower()
        if key not in allowed_sbf:
            raise ValueError(
                f"Invalid choice for sbf_v '{sbf_v}'. Must be one of: {allowed_sbf}"
            )

        sbf_v_column = sbf_map[key]
        if sbf_v_column not in sbf_hf_df.columns:
            raise KeyError(
                f"Column '{sbf_v_column}' not found in SBF Hubble-flow file. "
                f"Available columns include: {list(sbf_hf_df.columns)}"
            )



    methods, method_index, anchors, anchor_index, sources, source_index, hms, hms_index, mu_hms_error, mu_anchor_error = data_loader.compute_host_extras(host_df, anchor_names, anchor_g_mu_sigma)




    # Load other constraint files
    # Create paths only if the file is not "none"
    coma_path = os.path.join(datadir, coma_file) if coma_file else None
    mm_path = os.path.join(datadir, mm_file) if mm_file else None
    epm_path = os.path.join(datadir, epm_file) if epm_file else None
    coma_df = data_loader.load_coma(coma_path)
    mm_df = data_loader.load_mm(mm_path, q0, j0, vpec_error, vcorr_mm=vcorr_mm, mm_dz=mm_dz)# assuming distance-redshift case if mm_dz=True (default value)
    epm_df = data_loader.load_epm(epm_path, q0, j0, vpec_error, vcorr_epm=vcorr_epm, epm_dz=epm_dz)# assuming distance-redshift case if epm_dz=True (default value)


    # Build the system of equations
    eq_data = equations.build_equations(host_df, mm_df, epm_df, coma_df, groups_df,
                                          q0, j0, vpec_error, sbf_vpec_error,
                                          sn1a_calib_df=sn1a_calib_df, sn1a_hf_df=sn1a_hf_df, sn1a_hf_cov=sn1a_hf_cov, sn1a_vp_column=sn1a_vp_column, sn1a_ignore_offdiag=sn1a_ignore_offdiag, sn1a_IR=sn1a_IR,
                                          sn2_calib_df=sn2_calib_df, sn2_hf_df=sn2_hf_df, sn2_hf_cov=sn2_hf_cov,
                                          tf_calib_df=tf_calib_df, tf_hf_df=tf_hf_df, tf_hf_cov=tf_hf_cov, tf_v_column=tf_v_column,
                                          sbf_calib_df=sbf_calib_df, sbf_hf_df=sbf_hf_df, sbf_hf_cov=sbf_hf_cov, sbf_v_column=sbf_v_column,
                                          z_range_sn1a=z_range_sn1a, z_range_sn2=z_range_sn2, z_range_tf=z_range_tf, z_range_sbf=z_range_sbf,
                                          anchors=anchors, anchor_index=anchor_index, mu_anchor_error=mu_anchor_error,
                                          hms=hms, hms_index=hms_index, mu_hms_error=mu_hms_error,
                                          redshift_range=zbounds_hf,
                                          n_mas=n_mas, mas=mas, mas_index=mas_index, mas_error=mas_error,
                                          alpha_sn1a=alpha_sn1a, alpha_sn1a_error=alpha_sn1a_error,
                                          alpha_sn2=alpha_sn2, alpha_sn2_error=alpha_sn2_error,
                                          alpha_tf=alpha_tf, alpha_tf_error=alpha_tf_error,
                                          alpha_sbf=alpha_sbf, alpha_sbf_error=alpha_sbf_error,
                                          sigma_mag_fp=sigma_mag_fp)

    # Solve the system
    sol_results = solver.solve_system(eq_data)
    # Pass iabs and ihub for logging purposes
    sol_results["iabs"] = eq_data["iabs"]
    sol_results["ihub"] = eq_data["ihub"]


    # Log the results if requested
    if logfile and logfile.lower() != "none":
        logger.save_log(logfile, sol_results, host_exclude_list, calib_exclude)

 

    # Save comprehensive details to JSON if requested
    if savefile and savefile.lower() != "none":
        read_keys, skipped_keys = config_reader.give_keys()
        details = logger.generate_details(
            read_keys=read_keys,
            skipped_keys=skipped_keys,
            host_df=host_df,
            sn1a_calib_df=sn1a_calib_df,
            sn2_calib_df=sn2_calib_df,
            tf_calib_df=tf_calib_df,
            sbf_calib_df=sbf_calib_df,
            coma_df=coma_df,
            mm_df=mm_df,
            epm_df=epm_df,
            groups_df=groups_df,
            removed_hosts_without_calibrators=removed_hosts_without_calibrators,
            removed_sn1a_calib_without_hosts=removed_sn1a_calib_without_hosts,
            removed_sn2_calib_without_hosts=removed_sn2_calib_without_hosts,
            removed_tf_calib_without_hosts=removed_tf_calib_without_hosts,
            removed_sbf_calib_without_hosts=removed_sbf_calib_without_hosts,
            eq_data=eq_data,
            mas_ref_error=mas_ref_error,
            mas_ref_value=mas_ref_value,
            mas_ref_anchor=mas_ref_anchor,
            mas_ref_method=mas_ref_method,
            mas_ref_source=mas_ref_source,
            n_mas=n_mas,
            mas=mas,
            mas_index=mas_index,
            mas_error=mas_error,
            anchor_data=anchor_data,
            anchor_names=anchor_names,
            anchor_g_mu_sigma=anchor_g_mu_sigma,
            methods=methods,
            method_index=method_index,
            anchors=anchors,
            anchor_index=anchor_index,
            sources=sources,
            source_index=source_index,
            hms=hms,
            hms_index=hms_index,
            mu_hms_error=mu_hms_error,
            mu_anchor_error=mu_anchor_error,
            sol_results=sol_results
        )
        logger.save_details_as_json(savefile, details)

    # Print final results
    if verbose:
        print("\n=== Final Results ===")
        print(f" - log(H0) = {sol_results['logh0_value']:.5f} ± {np.sqrt(sol_results['logh0_var']):.5f}")
        print(f" - H0 = {sol_results['h0_value']:.5f} ± {sol_results['h0_error']:.5f} km/s/Mpc")
        if sol_results["iabs"] is not None:
            print(f" - MB = {sol_results['mzero_value']:.3f} ± {sol_results['mzero_error']:.3f}")
        if sol_results["ndof"] > 0:
            print(f" - Chi-squared: {sol_results['chi2']:.4f}, Degrees of Freedom: {sol_results['ndof']}, Number of Parameters: {sol_results['npars']}")
            print(f" - Reduced Chi-squared: {sol_results['chi2']/sol_results['ndof']:.4f}")
        else:
            print(" - No degrees of freedom.")
    else:
        print(f"H0={sol_results['h0_value']:.4f}+/-{sol_results['h0_error']:.4f} km/s/Mpc | Chi-squared={sol_results['chi2']:.4f} | ndof={sol_results['ndof']} | Reduced Chi-squared={('{:.4f}'.format(sol_results['chi2']/sol_results['ndof']) if sol_results['ndof']>0 else 'n/a')} | MB={(('{:.3f}+/-{:.3f}'.format(sol_results['mzero_value'], sol_results['mzero_error'])) if sol_results.get('iabs') is not None else 'n/a')}")




if __name__ == "__main__":
    main()
