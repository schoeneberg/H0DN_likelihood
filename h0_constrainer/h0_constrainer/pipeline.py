import os
import sys
import numpy as np
import pandas as pd
import scipy.linalg
from h0_constrainer import config_reader, data_loader, equations, solver, logger

class H0DN_Pipeline:
  ini = None
  data = None
  eq_data = None
  sol_results = None
  
  test = None
  
  def read_config(self, config_file):
    config_reader.reset_keys()
    
    ini = {}
    config = config_reader.load_config(config_file)
    config_dir = os.path.dirname(os.path.abspath(config_file))
    raw_datadir = config_reader.get_config_value(config, "DEFAULT", "datadir", "./")
    ini['datadir'] = os.path.abspath(os.path.join(config_dir, raw_datadir)) # IGNORES config_dir if raw_datadir is an aboslute path, see documentation of os.path.join
    ini['verbose'] = config_reader.get_config_value(config, "DEFAULT", "verbose", False, as_bool=True)
    ini['warn'] = config_reader.get_config_value(config, "DEFAULT", "warn", True, as_bool=True)

    config_reader.set_flags(verbose_flag=ini['verbose'], warn_flag=ini['warn'])

    ini['anchor_file'] = config_reader.get_config_value(config, "DEFAULT", "anchors", "anchors.dat").strip().lower()
    ini['host_file'] = config_reader.get_config_value(config, "DEFAULT", "host_data", "").strip().lower()

    ini['sn1a_calib_file'] = config_reader.get_config_value(config, "DEFAULT", "sn1a_calib", "").strip().lower()
    ini['sn2_calib_file'] = config_reader.get_config_value(config, "DEFAULT", "sn2_calib", "").strip().lower()
    ini['tf_calib_file'] = config_reader.get_config_value(config, "DEFAULT", "tf_calib", "").strip().lower()
    ini['sbf_calib_file'] = config_reader.get_config_value(config, "DEFAULT", "sbf_calib", "").strip().lower()

    ini['groups_file'] = config_reader.get_config_value(config, "DEFAULT", "groups", "").strip().lower()

    ini['sn1a_hf_file'] = config_reader.get_config_value(config, "DEFAULT", "sn1a_hf", "").strip().lower()
    ini['sn1a_hf_cov_file'] = config_reader.get_config_value(config, "DEFAULT", "sn1a_hf_cov", "").strip().lower()
    ini['sn2_hf_file'] = config_reader.get_config_value(config, "DEFAULT", "sn2_hf", "").strip().lower()
    ini['sn2_hf_cov_file'] = config_reader.get_config_value(config, "DEFAULT", "sn2_hf_cov", "").strip().lower()
    ini['sn1a_ignore_offdiag'] = config_reader.get_config_value(config, "DEFAULT", "sn1a_ignore_offdiag", False, as_bool=True)
    ini['sn1a_IR'] = config_reader.get_config_value(config, "DEFAULT", "sn1a_IR", False, as_bool=True)
    ini['sn1a_intrinsic'] = config_reader.get_config_value(config, "DEFAULT", "sn1a_intrinsic", 0.0, as_float=True)
    ini['tf_hf_file'] = config_reader.get_config_value(config, "DEFAULT", "tf_hf", "").strip().lower()
    ini['tf_hf_cov_file'] = config_reader.get_config_value(config, "DEFAULT", "tf_hf_cov", "").strip().lower()
    ini['sbf_hf_file'] = config_reader.get_config_value(config, "DEFAULT", "sbf_hf", "").strip().lower()
    ini['sbf_hf_cov_file'] = config_reader.get_config_value(config, "DEFAULT", "sbf_hf_cov", "").strip().lower()

    ini['mm_file'] = config_reader.get_config_value(config, "DEFAULT", "mm", "").strip().lower()
    ini['epm_file'] = config_reader.get_config_value(config, "DEFAULT", "epm", "").strip().lower()    
    ini['coma_file'] = config_reader.get_config_value(config, "DEFAULT", "coma", "").strip().lower()
    ini['sigma_mag_fp'] = config_reader.get_config_value(config, "DEFAULT", "sigma_mag_fp", 0.0371, as_float=True)

    ini['mm_dz'] = config_reader.get_config_value(config, "DEFAULT", "mm_dz", True, as_bool=True)
    ini['epm_dz'] = config_reader.get_config_value(config, "DEFAULT", "epm_dz", True, as_bool=True)

    ini['logfile'] = config_reader.get_config_value(config, "DEFAULT", "logfile", "none").strip()
    ini['savefile'] = config_reader.get_config_value(config, "DEFAULT", "savefile", "none").strip()
    # Exclusion lists
    host_exclude_list_str = config_reader.get_config_value(config, "DEFAULT", "host_exclude_list", "")
    ini['host_exclude_list'] = [item.strip() for item in host_exclude_list_str.split(",")] if host_exclude_list_str else []
    calib_exclude_str = config_reader.get_config_value(config, "DEFAULT", "calib_exclude", "")
    ini['calib_exclude'] = [item.strip() for item in calib_exclude_str.split(",")] if calib_exclude_str else []
    # Inclusion lists
    host_include_list_str = config_reader.get_config_value(config, "DEFAULT", "host_include_list", "")
    ini['host_include_list'] = [item.strip() for item in host_include_list_str.split(",")] if host_include_list_str else []
    calib_include_str = config_reader.get_config_value(config, "DEFAULT", "calib_include", "")
    ini['calib_include'] = [item.strip() for item in calib_include_str.split(",")] if calib_include_str else []

    ini['q0'] = config_reader.get_config_value(config, "DEFAULT", "q0", -0.55, as_float=True)
    ini['j0'] = config_reader.get_config_value(config, "DEFAULT", "j0", 1.0, as_float=True)

    zmin_hf = config_reader.get_config_value(config, "DEFAULT", "zmin_hf", 0.0, as_float=True)
    zmax_hf = config_reader.get_config_value(config, "DEFAULT", "zmax_hf", 100.0, as_float=True)
    if zmin_hf == 0.0 and zmax_hf == 100.0:
        zbounds_hf = None
    else:
        zbounds_hf = (zmin_hf, zmax_hf)
    ini['zbounds_hf'] = zbounds_hf
    ini['z_range_sn1a'] = config_reader.get_config_value(config, "DEFAULT", "z_range_sn1a", False, as_bool=True)
    ini['z_range_sn2'] = config_reader.get_config_value(config, "DEFAULT", "z_range_sn2", False, as_bool=True)
    ini['z_range_tf'] = config_reader.get_config_value(config, "DEFAULT", "z_range_tf", False, as_bool=True)
    ini['z_range_sbf'] = config_reader.get_config_value(config, "DEFAULT", "z_range_sbf", False, as_bool=True)

    ini['alpha_sn1a'] = config_reader.get_config_value(config, "DEFAULT", "alpha_sn1a", "", as_float=True)
    ini['alpha_sn1a_error'] = config_reader.get_config_value(config, "DEFAULT", "alpha_sn1a_error", "", as_float=True)
    ini['alpha_sn2'] = config_reader.get_config_value(config, "DEFAULT", "alpha_sn2", "", as_float=True)
    ini['alpha_sn2_error'] = config_reader.get_config_value(config, "DEFAULT", "alpha_sn2_error", "", as_float=True)
    ini['alpha_tf'] = config_reader.get_config_value(config, "DEFAULT", "alpha_tf", "", as_float=True)
    ini['alpha_tf_error'] = config_reader.get_config_value(config, "DEFAULT", "alpha_tf_error", "", as_float=True)
    ini['alpha_sbf'] = config_reader.get_config_value(config, "DEFAULT", "alpha_sbf", "", as_float=True)
    ini['alpha_sbf_error'] = config_reader.get_config_value(config, "DEFAULT", "alpha_sbf_error", "", as_float=True)

    ini['unused_hosts_via_covariance'] = config_reader.get_config_value(config, "DEFAULT", "unused_hosts_via_covariance", False, as_bool=True)

    ini['vcorr_mm'] = config_reader.get_config_value(config, "DEFAULT", "vcorr_mm", "2M++")
    ini['vcorr_epm'] = config_reader.get_config_value(config, "DEFAULT", "vcorr_epm", "2M++")
    ini['vpec_error'] = config_reader.get_config_value(config, "DEFAULT", "vpec_error", 240.0, as_float=True)
    ini['sbf_vpec_error'] = config_reader.get_config_value(config, "DEFAULT", "sbf_vpec_error", ini['vpec_error'], as_float=True)
    ini['sn1a_vp'] = config_reader.get_config_value(config, "DEFAULT", "sn1a_vp", "2m++")
    ini['tf_v'] = config_reader.get_config_value(config, "DEFAULT", "tf_v", "2m++")
    ini['sbf_v'] = config_reader.get_config_value(config, "DEFAULT", "sbf_v", "2m++")

    if ini['verbose']:
        config_reader.print_config_keys()

    read, skipped = config_reader.give_keys()
    read = [r[0] for r in read]
  
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
    
    self.ini = ini
    return ini

  def read_data(self, ini = None):
    ini = self.check_for(ini, "ini", "read_data")
    
    data = {}
    # TODO :: REMOVE
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

    datadir = ini['datadir']
    # Load data files
    host_df = data_loader.load_hosts(datadir, ini['host_file'])
    if ini['sn1a_IR']:
        sn1a_calib_df = data_loader.load_calibrators(datadir, ini['sn1a_calib_file'], "sn_IR")
    else:
        sn1a_calib_df = data_loader.load_calibrators(datadir, ini['sn1a_calib_file'], "sn")
    if sn1a_calib_df is not None:
        sn1a_calib_df["sigma"] = np.sqrt(sn1a_calib_df["sigma"]**2 + ini['sn1a_intrinsic']**2)

    sn2_calib_df = data_loader.load_calibrators(datadir, ini['sn2_calib_file'], "sn")
    tf_calib_df = data_loader.load_calibrators(datadir, ini['tf_calib_file'], "tf")
    sbf_calib_df = data_loader.load_calibrators(datadir, ini['sbf_calib_file'], "sbf")

    groups_df = data_loader.load_groups(datadir, ini['groups_file'])

    # Host data filtering (apply inclusion/exclusion lists)
    if host_df is not None:
        host_include_list, host_exclude_list = ini['host_include_list'], ini['host_exclude_list']
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

    calib_include, calib_exclude = ini['calib_include'], ini['calib_exclude']

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
    if ini['verbose']:
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
    if not ini['unused_hosts_via_covariance']:
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
    anchor_data, anchor_names, anchor_g_mu_sigma = data_loader.load_anchors(datadir, anchorfile=ini['anchor_file'])


    # --- SN Ia load Hubble flow ---
    sn1a_vp_column = None
    if ini['sn1a_hf_file']:
        if ini['sn1a_hf_cov_file']:
            sn1a_hf_df, sn1a_hf_cov = data_loader.load_hf(
                datadir, ini['sn1a_hf_file'], "sn1a", covar_filename=ini['sn1a_hf_cov_file']
            )
        else:
            method = "sn1a_IR" if ini['sn1a_IR'] else "sn1a"
            sn1a_hf_df, sn1a_hf_cov = data_loader.load_hf(datadir, ini['sn1a_hf_file'], method)
        sn1a_hf_df["mb_err"] = np.sqrt(sn1a_hf_df["mb_err"]**2 + ini['sn1a_intrinsic']**2)

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

        key = ini['sn1a_vp'].strip().lower()

        if not ini['sn1a_IR']:
            # ---- Non-IR: use full SN Ia options (including 'plain') ----
            if key not in allowed_sn1a_usual:
                raise ValueError(
                    f"Invalid choice for sn1a_vp '{ini['sn1a_vp']}'. Must be one of: {allowed_sn1a_usual}"
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
    if ini['sn2_hf_file']:
        sn2_hf_df, sn2_hf_cov = data_loader.load_hf(
            datadir, ini['sn2_hf_file'], "sn2", covar_filename=ini['sn2_hf_cov_file']
        )

    # --- Tully–Fisher load Hubble flow ---
    tf_v_column = None
    if ini['tf_hf_file']:
        tf_hf_df, tf_hf_cov = data_loader.load_hf(
            datadir, ini['tf_hf_file'], "tf", covar_filename=ini['tf_hf_cov_file']
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
    if ini['sbf_hf_file']:
        sbf_hf_df, sbf_hf_cov = data_loader.load_hf(
            datadir, ini['sbf_hf_file'], "sbf", covar_filename=ini['sbf_hf_cov_file']
        )

        # SBF supports only: 'cmb' (mapped to group frame) and '2m++'
        allowed_sbf = ["cmb", "2m++"]
        sbf_map = {
            "cmb":  "vgrp",   # 'cmb' input -> group-frame column
            "2m++": "v2m++",
        }

        key = ini['sbf_v'].strip().lower()
        if key not in allowed_sbf:
            raise ValueError(
                f"Invalid choice for sbf_v '{ini['sbf_v']}'. Must be one of: {allowed_sbf}"
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
    coma_path = os.path.join(datadir, ini['coma_file']) if ini['coma_file'] else None
    mm_path = os.path.join(datadir, ini['mm_file']) if ini['mm_file'] else None
    epm_path = os.path.join(datadir, ini['epm_file']) if ini['epm_file'] else None
    coma_df = data_loader.load_coma(coma_path)
    mm_df = data_loader.load_mm(mm_path, ini['q0'], ini['j0'], ini['vpec_error'], vcorr_mm=ini['vcorr_mm'], mm_dz=ini['mm_dz'])# assuming distance-redshift case if mm_dz=True (default value)
    epm_df = data_loader.load_epm(epm_path, ini['q0'], ini['j0'], ini['vpec_error'], vcorr_epm=ini['vcorr_epm'], epm_dz=ini['epm_dz'])# assuming distance-redshift case if epm_dz=True (default value)
    
    data = {'host':host_df,'mm':mm_df,'epm':epm_df,'coma':coma_df,'groups':groups_df,
            'sn1a_calib':sn1a_calib_df,'sn1a_hf':sn1a_hf_df,
            'sn2_calib':sn2_calib_df,'sn2_hf':sn2_hf_df,
            'tf_calib':tf_calib_df,'tf_hf':tf_hf_df,
            'sbf_calib':sbf_calib_df,'sbf_hf':sbf_hf_df,
            'anchors':anchors,'anchor_index':anchor_index,'mu_anchor_error':mu_anchor_error,
            'hms':hms,'hms_index':hms_index,'mu_hms_error':mu_hms_error,
            'n_mas':n_mas,'mas':mas,'mas_index':mas_index, 'mas_error':mas_error,
            'sn1a_hf_cov':sn1a_hf_cov, 'sn1a_vp_column':sn1a_vp_column,
            'sn2_hf_cov':sn2_hf_cov,
            'tf_hf_cov':tf_hf_cov,'tf_v_column':tf_v_column,
            'sbf_hf_cov':sbf_hf_cov,'sbf_v_column':sbf_v_column,
            'removed_hosts_without_calibrators':removed_hosts_without_calibrators, 'removed_sn1a_calib_without_hosts':removed_sn1a_calib_without_hosts, 'removed_sn2_calib_without_hosts':removed_sn2_calib_without_hosts, 'removed_tf_calib_without_hosts':removed_tf_calib_without_hosts, 'removed_sbf_calib_without_hosts':removed_sbf_calib_without_hosts,
            'mas_ref_error':mas_ref_error, 'mas_ref_value':mas_ref_value, 'mas_ref_anchor':mas_ref_anchor, 'mas_ref_method':mas_ref_method, 'mas_ref_source':mas_ref_source,
            'anchor_data':anchor_data, 'anchor_names':anchor_names, 'anchor_g_mu_sigma':anchor_g_mu_sigma,
            'methods':methods,'method_index':method_index,'sources':sources,'source_index':source_index}

    self.data = data
    return data
    
  def build_equations(self, ini = None, data = None):
    ini = self.check_for(ini, "ini", "build_equations")
    data = self.check_for(data, "data", "build_equations")

    eq_data = equations.build_equations(data['host'], data['mm'], data['epm'], data['coma'], data['groups'],
                                        ini['q0'], ini['j0'], ini['vpec_error'], ini['sbf_vpec_error'],
                                        sn1a_calib_df=data['sn1a_calib'], sn1a_hf_df=data['sn1a_hf'], sn1a_hf_cov=data['sn1a_hf_cov'], sn1a_vp_column=data['sn1a_vp_column'], sn1a_ignore_offdiag=ini['sn1a_ignore_offdiag'], sn1a_IR=ini['sn1a_IR'],
                                        sn2_calib_df=data['sn2_calib'], sn2_hf_df=data['sn2_hf'], sn2_hf_cov=data['sn2_hf_cov'],
                                        tf_calib_df=data['tf_calib'], tf_hf_df=data['tf_hf'], tf_hf_cov=data['tf_hf_cov'], tf_v_column=data['tf_v_column'],
                                        sbf_calib_df=data['sbf_calib'], sbf_hf_df=data['sbf_hf'], sbf_hf_cov=data['sbf_hf_cov'], sbf_v_column=data['sbf_v_column'],
                                        z_range_sn1a=ini['z_range_sn1a'], z_range_sn2=ini['z_range_sn2'], z_range_tf=ini['z_range_tf'], z_range_sbf=ini['z_range_sbf'],
                                        anchors=data['anchors'], anchor_index=data['anchor_index'], mu_anchor_error=data['mu_anchor_error'],
                                        hms=data['hms'], hms_index=data['hms_index'], mu_hms_error=data['mu_hms_error'],
                                        redshift_range=ini['zbounds_hf'],
                                        n_mas=data['n_mas'], mas=data['mas'], mas_index=data['mas_index'], mas_error=data['mas_error'],
                                        alpha_sn1a=ini['alpha_sn1a'], alpha_sn1a_error=ini['alpha_sn1a_error'],
                                        alpha_sn2=ini['alpha_sn2'], alpha_sn2_error=ini['alpha_sn2_error'],
                                        alpha_tf=ini['alpha_tf'], alpha_tf_error=ini['alpha_tf_error'],
                                        alpha_sbf=ini['alpha_sbf'], alpha_sbf_error=ini['alpha_sbf_error'],
                                        sigma_mag_fp=ini['sigma_mag_fp'])
    self.eq_data = eq_data
    return eq_data

  def adjust_equations(self, cosmo, ini = None, data = None, eq_data = None):
    eq_data = self.check_for(eq_data, "eq_data", "adjust_equations")
    data = self.check_for(data, "data", "adjust_equations")
    ini = self.check_for(ini, "ini", "adjust_equations")
    eq_data = equations.adjust_equation_data(cosmo, eq_data, data['mm'], data['epm'], data['sbf_hf'], data['tf_hf'], data['sn2_hf'], data['sn1a_hf'],
                         redshift_range = ini['zbounds_hf'], vpec_error = ini['vpec_error'], sbf_vpec_error = ini['sbf_vpec_error'],
                         sn1a_IR = ini['sn1a_IR'], sn1a_vp_column = data['sn1a_vp_column'], sn1a_hf_cov = data['sn1a_hf_cov'], z_range_sn1a = ini['z_range_sn1a'], sn1a_ignore_offdiag = ini['sn1a_ignore_offdiag'],
                         z_range_sn2 = ini['z_range_sn2'],
                         tf_v_column = data['tf_v_column'], z_range_tf = ini['z_range_tf'],
                         sbf_v_column = data['sbf_v_column'], z_range_sbf = ini['z_range_sbf'],
                         alpha_sn1a=ini['alpha_sn1a'], alpha_sn1a_error=ini['alpha_sn1a_error'],
                         alpha_sn2=ini['alpha_sn2'], alpha_sn2_error=ini['alpha_sn2_error'],
                         alpha_tf=ini['alpha_tf'], alpha_tf_error=ini['alpha_tf_error'],
                         alpha_sbf=ini['alpha_sbf'], alpha_sbf_error=ini['alpha_sbf_error'])
    self.eq_data = eq_data
    return eq_data

  def solve(self, eq_data = None):
    eq_data = self.check_for(eq_data, "eq_data", "solve")

    # Solve the system
    sol_results = solver.solve_system(eq_data)
    # Pass iabs and ihub for logging purposes
    sol_results["iabs"] = eq_data["iabs"]
    sol_results["ihub"] = eq_data["ihub"]

    self.sol_results = sol_results
    return sol_results

  def log_results(self, ini = None, sol_results = None):
    sol_results = self.check_for(sol_results, "sol_results", "log_results")
    ini = self.check_for(ini, "ini", "log_results")

    if ini['logfile'] and ini['logfile'].lower() != "none":
      logger.save_log(ini['logfile'], sol_results, ini['host_exclude_list'], ini['calib_exclude'])

  def save_results(self, ini = None, data = None, eq_data = None, sol_results = None):
    sol_results = self.check_for(sol_results, "sol_results", "save_results")
    eq_data = self.check_for(eq_data, "eq_data", "save_results")
    data = self.check_for(data, "data", "save_results")
    ini = self.check_for(ini, "ini", "save_results")

    # Save comprehensive details to JSON if requested
    if ini['savefile'] and ini['savefile'].lower() != "none":
        read_keys, skipped_keys = config_reader.give_keys()
        details = logger.generate_details(
            read_keys=read_keys,
            skipped_keys=skipped_keys,
            host_df=data['host'],
            sn1a_calib_df=data['sn1a_calib'],
            sn2_calib_df=data['sn2_calib'],
            tf_calib_df=data['tf_calib'],
            sbf_calib_df=data['sbf_calib'],
            coma_df=data['coma'],
            mm_df=data['mm'],
            epm_df=data['epm'],
            groups_df=data['groups'],
            removed_hosts_without_calibrators=data['removed_hosts_without_calibrators'],
            removed_sn1a_calib_without_hosts=data['removed_sn1a_calib_without_hosts'],
            removed_sn2_calib_without_hosts=data['removed_sn2_calib_without_hosts'],
            removed_tf_calib_without_hosts=data['removed_tf_calib_without_hosts'],
            removed_sbf_calib_without_hosts=data['removed_sbf_calib_without_hosts'],
            eq_data=eq_data,
            mas_ref_error=data['mas_ref_error'],
            mas_ref_value=data['mas_ref_value'],
            mas_ref_anchor=data['mas_ref_anchor'],
            mas_ref_method=data['mas_ref_method'],
            mas_ref_source=data['mas_ref_source'],
            n_mas=data['n_mas'],
            mas=data['mas'],
            mas_index=data['mas_index'],
            mas_error=data['mas_error'],
            anchor_data=data['anchor_data'],
            anchor_names=data['anchor_names'],
            anchor_g_mu_sigma=data['anchor_g_mu_sigma'],
            methods=data['methods'],
            method_index=data['method_index'],
            anchors=data['anchors'],
            anchor_index=data['anchor_index'],
            sources=data['sources'],
            source_index=data['source_index'],
            hms=data['hms'],
            hms_index=data['hms_index'],
            mu_hms_error=data['mu_hms_error'],
            mu_anchor_error=data['mu_anchor_error'],
            sol_results=sol_results
        )
        logger.save_details_as_json(ini['savefile'], details)
  
  def print_results(self, ini = None, sol_results = None):
    ini = self.check_for(ini, "ini", "print_results")
    sol_results = self.check_for(sol_results, "sol_results", "print_results")
  
    # Print final results
    if ini['verbose']:
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





  # Note that since the likelihood is Gaussian, the option
  # 'minimize_except' is the same as 'marginalize_except',
  # up to a normalization constant
  def get_likelihood(self, eq_data = None, inv_covar = None, minimize_except="all", verbose=False):
    eq_data = self.check_for(eq_data, "eq_data", "get_likelihood")
    if inv_covar is None:
      # About as fast as inverting the covar from the eq_data manually.
      # Should in any case only be done rarely, as the inv_covar should be provided if speed is of interest
      sol_results = self.solve(eq_data = eq_data) 
      inv_covar = sol_results['inv_covar']

    # Let's compute the likelihood
    # We simplify it to const + dX^T  gamma  dX
    # This is equivalent to the normal form
    # (Y - A*X)^T C^(-1)  (Y-A*X)
    # but a whole lot faster
    # Basically this works by factoring out the X and then completing the square

    inv_C = inv_covar
    Y = eq_data["yval"]
    A = eq_data["coeffs"]
    pars = [x['objname'] for x in eq_data["parinfo"]]

    # Now, let's bake
    alpha = Y.T @ inv_C @ Y
    beta = Y.T @ inv_C @ A
    gamma = A.T @ inv_C @ A

    # Get the analytical bestfit
    X_bf = np.linalg.solve(gamma, beta)
    chi2_min = alpha - beta @ X_bf #+ (eq_data["chisq_sn1a_hf"] if eq_data["chisq_sn1a_hf"] else 0)
    # TODO :: Possibly add the chi^2 contribtions from the computation of the alpha intercepts, if asked !
    X_cov = np.linalg.inv(gamma)

    if minimize_except!="all":
      mask = np.array([any([(p in x) for p in minimize_except]) for x in pars],dtype=bool)
      pars = [pars[i] for i in range(len(pars)) if mask[i]] #doing it the old fashioned way, since pars is a list
      X_bf = X_bf[mask]
      g11 = gamma[mask][:, mask]
      g12 = gamma[mask][:, ~mask]
      g22 = gamma[~mask][:, ~mask]
      gamma = g11 - g12 @ np.linalg.inv(g22) @ g12.T
      X_cov = np.linalg.inv(gamma)

    X_bf_vec = X_bf[np.newaxis, :]

    def log_likelihood(paramvector):
      """
      For a parameter vector of shape (k,p) returns (k,) vector of chi^2 values.
      """
      # Ensure paramvector is 2d but in a way that will not slow down the code at all
      paramvector = paramvector if paramvector.ndim > 1 else paramvector[np.newaxis, :]
      dX = paramvector - X_bf_vec
      # The last part is some super fast numpy magic to do X^T @ gamma @ X in parallel
      return -0.5*(chi2_min + np.sum(dX @ gamma * dX,axis=1))
    return pars, log_likelihood, X_bf, X_cov

  def check_for(self, new_state, name_of_state, name_of_function):
    # These are the functions that are required to compute a given internal state of name "name_of_state"
    function = {"ini":"read_config", "data":"read_data","eq_data":"build_equations", "sol_results":"solve"}
    if new_state is None:
      if getattr(self,name_of_state) is None:
        raise ValueError(f"Need to call '{function[name_of_state]}' before '{name_of_function}' when not explicitly passing '{name_of_state}'.")
      else:
        new_state = getattr(self, name_of_state)
    return new_state


