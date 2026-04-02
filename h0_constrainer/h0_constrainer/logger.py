"""Logging and results serialization module.

This module handles saving analysis results to log files and JSON format.
It provides utilities for converting complex data structures (including numpy
arrays and pandas DataFrames) into JSON-serializable formats.

Functions:
    save_log: Save basic results summary to a text log file.
    generate_details: Create comprehensive results dictionary with all analysis details.
    save_details_as_json: Save detailed results to JSON format.
    make_json_serializable: Convert numpy/pandas objects to JSON-compatible types.
"""

import json
import os
import numpy as np
import pandas as pd
from h0_constrainer import config_reader

def save_log(logfile, results, exclude_list, calib_exclude):
    """Save basic analysis results to a text log file.
    
    Writes a summary including H0 value, absolute magnitude (if applicable),
    chi-squared statistics, and lists of excluded objects.
    
    Args:
        logfile (str): Path to output log file.
        results (dict): Results dictionary from solve_system() containing:
            - h0_value, h0_error: Hubble constant and uncertainty.
            - iabs: Index of absolute magnitude parameter (or None).
            - params: Array of fitted parameters.
            - invsolmat: Inverse solution matrix for uncertainties.
            - chi2: Chi-squared value.
            - ndof: Degrees of freedom.
        exclude_list (list): List of excluded hosts.
        calib_exclude (list): List of excluded calibrators.
    """
    with open(logfile, "w") as log:
        log.write("=== RESULTS ===\n")
        log.write(f"H0 Value and uncertainty: {results['h0_value']:.5f} ± {results['h0_error']:.5f}\n")
        if results["iabs"] is not None:
            mzero_value = results["params"][results["iabs"]]
            mzero_error = (results["invsolmat"][results["iabs"], results["iabs"]]) ** 0.5
            log.write(f"M_B value and uncertainty: {mzero_value:.5f} ± {mzero_error:.5f}\n")
        log.write(f"Chi²: {results['chi2']:.5f}, Degrees of Freedom: {results['ndof']}\n")
        log.write(f"Excluded Hosts: {exclude_list}\n")
        log.write(f"Excluded Calibrators: {calib_exclude}\n")
    config_reader.vprint(f"\nResults successfully saved to log file: {logfile}")



def generate_details(
    read_keys,
    skipped_keys,
    # --- Data frames ---
    host_df,
    sn1a_calib_df,
    sn2_calib_df,
    tf_calib_df,
    sbf_calib_df,
    coma_df,
    mm_df,
    epm_df,
    groups_df,
    # --- From rung_equalize ---
    removed_hosts_without_calibrators,
    removed_sn1a_calib_without_hosts,
    removed_sn2_calib_without_hosts,
    removed_tf_calib_without_hosts,
    removed_sbf_calib_without_hosts,
    # --- From build_equations ---
    eq_data,
    # --- From drop_anchor_rows ---
    mas_ref_error,
    mas_ref_value,
    mas_ref_anchor,
    mas_ref_method,
    mas_ref_source,
    n_mas,
    mas,
    mas_index,
    mas_error,
    # --- From load_anchors ---
    anchor_data,
    anchor_names,
    anchor_g_mu_sigma,
    # --- From compute_host_extras ---
    methods,
    method_index,
    anchors,
    anchor_index,
    sources,
    source_index,
    hms,
    hms_index,
    mu_hms_error,
    mu_anchor_error,
    # --- From solve_system ---
    sol_results
):
    """Generate comprehensive details dictionary for JSON serialization.

    Aggregates all analysis information including configuration, data, equation system,
    and results into a single nested dictionary structure suitable for saving to JSON.

    Args:
        read_keys (list): Configuration keys that were read.
        skipped_keys (list): Configuration keys that were skipped.
        host_df (pd.DataFrame): Host distance measurements.
        sn1a_calib_df (pd.DataFrame): SNe Ia calibrator data.
        sn2_calib_df (pd.DataFrame): SNe II calibrator data.
        tf_calib_df (pd.DataFrame): Tully-Fisher calibrator data.
        sbf_calib_df (pd.DataFrame): SBF calibrator data.
        coma_df (pd.DataFrame): Coma cluster data.
        mm_df (pd.DataFrame): Megamaser data.
        epm_df (pd.DataFrame): EPM SN II data.
        groups_df (pd.DataFrame): Galaxy group data.
        eq_data (dict): Equation system data from build_equations().
        mas_ref_error, mas_ref_value, mas_ref_anchor, mas_ref_method, mas_ref_source:
            Method-Anchor-Source (MAS) reference arrays from drop_anchor_rows().
        n_mas (int): Number of unique MAS combinations.
        mas (np.ndarray): Unique MAS names.
        mas_index (np.ndarray): MAS index for each host row.
        mas_error (np.ndarray): MAS uncertainties.
        anchor_data (pd.DataFrame): Anchor distance data from load_anchors().
        anchor_names (np.ndarray): Anchor names.
        anchor_g_mu_sigma (np.ndarray): Anchor distance modulus uncertainties.
        methods, method_index, anchors, anchor_index, sources, source_index:
            Unique methods, anchors, sources and their index mappings from compute_host_extras().
        hms, hms_index, mu_hms_error: Host-Method-Source combinations and uncertainties.
        mu_anchor_error (np.ndarray): Anchor uncertainties.
        sol_results (dict): Solution results from solve_system().

    Returns:
        dict: Nested dictionary containing all analysis details organized by category.
    """
    read_dict = dict(read_keys)

    details = {
        # =====================================================================
        # CONFIGURATION
        # =====================================================================
        "config_info": {
            "read_keys": read_keys,
            "skipped_keys": skipped_keys
        },

        # =====================================================================
        # PRIMARY RESULTS
        # =====================================================================
        "h0_value": sol_results.get("h0_value"),
        "h0_error": sol_results.get("h0_error"),
        "logh0_value": sol_results.get("logh0_value"),
        "logh0_var": sol_results.get("logh0_var"),
        "chi2": sol_results.get("chi2"),
        "ndof": sol_results.get("ndof"),
        "ndof_full": sol_results.get("ndof_full"),
        "npars": sol_results.get("npars"),
        "mzero_value": sol_results.get("mzero_value"),
        "mzero_error": sol_results.get("mzero_error"),
        "mu_coma_value": sol_results.get("mu_coma_value"),
        "mu_coma_error": sol_results.get("mu_coma_error"),

        # =====================================================================
        # EQUATION SYSTEM STRUCTURE
        # =====================================================================
        "parinfo": eq_data.get("parinfo"),
        "eq_descr": eq_data.get("eq_descr"),
        "eq_shape": eq_data.get("eq_shape"),
        "eq_covar": eq_data.get("eq_covar"),
        "coeffs": eq_data.get("coeffs"),
        "yval": eq_data.get("yval"),
        "residuals": sol_results.get("residuals"),

        # Parameter indices
        "ihub": sol_results.get("ihub"),
        "iabs": sol_results.get("iabs"),
        "iabs_sn2": eq_data.get("iabs_sn2"),
        "iabs_tf": eq_data.get("iabs_tf"),
        "iabs_sbf": eq_data.get("iabs_sbf"),
        "icoma": eq_data.get("icoma"),

        # Equation start indices
        "ieq_sn1a_start": eq_data.get("ieq_sn1a_start"),
        "ieq_sn2_start": eq_data.get("ieq_sn2_start"),
        "ieq_tf_start": eq_data.get("ieq_tf_start"),
        "ieq_sbf_start": eq_data.get("ieq_sbf_start"),
        "ieq_coma_start": eq_data.get("ieq_coma_start"),
        "ieq_mm_start": eq_data.get("ieq_mm_start"),
        "ieq_epm_start": eq_data.get("ieq_epm_start"),
        "ieq_h0_m1a": eq_data.get("ieq_h0_m1a"),
        "ieq_h0_m2": eq_data.get("ieq_h0_m2"),
        "ieq_h0_mtf": eq_data.get("ieq_h0_mtf"),
        "ieq_h0_msbf": eq_data.get("ieq_h0_msbf"),

        # =====================================================================
        # ANCHOR DETAILS
        # =====================================================================
        "anchor_details": {
            "n_anchors": len(anchors),
            "anchors": anchors,
            "mu_anchor_error": mu_anchor_error,
            "anchor_names": anchor_names,
            "data_mu_value": anchor_data,
            "anchor_g_mu_sigma": anchor_g_mu_sigma
        },

        # =====================================================================
        # MAS (METHOD-ANCHOR-SOURCE) DETAILS
        # =====================================================================
        "mas_details": {
            "n_mas": n_mas,
            "mas": mas,
            "mas_error": mas_error,
            "mas_ref_name": mas,
            "mas_ref_value": mas_ref_value,
            "mas_ref_error": mas_ref_error,
            "mas_ref_anchor": mas_ref_anchor,
            "mas_ref_method": mas_ref_method,
            "mas_ref_source": mas_ref_source
        },

        # =====================================================================
        # HOST DATA DETAILS
        # =====================================================================
        "host_details": {
            "n_hosts": len(host_df),
            "hosts": host_df["host"],
            "methodname": [methods[i] for i in method_index],
            "anchorname": [anchors[i] for i in anchor_index],
            "sourcename": [sources[i] for i in source_index],
            "hms_name": hms,
            "hms_index": hms_index,
            "hms_error": mu_hms_error,
            "mu_host_value": host_df["mu_host"],
            "mu_host_error": host_df["mu_error"],
            "mas_name": mas,
            "mas_index": mas_index,
            "anchor_host_error": eq_data.get("anchor_host_error"),
            "mas_host_error": eq_data.get("mas_host_error")
        },

        # =====================================================================
        # SNe Ia CALIBRATOR DETAILS
        # =====================================================================
        "sn1a_calib_details": {
            "n_sn1a": len(sn1a_calib_df) if sn1a_calib_df is not None else 0,
            "hosts": sn1a_calib_df["host"] if sn1a_calib_df is not None else [],
            "sn_names": sn1a_calib_df["SN"] if sn1a_calib_df is not None else [],
            "magnitudes": sn1a_calib_df["m0_Bi"] if sn1a_calib_df is not None else [],
            "sigma": sn1a_calib_df["sigma"] if sn1a_calib_df is not None else [],
            # Alpha from Hubble flow fit
            "a_sn1a": eq_data.get("a_sn1a"),
            "a_sn1a_err": eq_data.get("a_sn1a_err"),
            "chisq_sn1a_hf": eq_data.get("chisq_sn1a_hf"),
            "ndof_sn1a_hf": eq_data.get("ndof_sn1a_hf"),
        } if sn1a_calib_df is not None and len(sn1a_calib_df) > 0 else {},

        # =====================================================================
        # SNe II CALIBRATOR DETAILS
        # =====================================================================
        "sn2_calib_details": {
            "n_sn2": len(sn2_calib_df) if sn2_calib_df is not None else 0,
            "hosts": sn2_calib_df["host"] if sn2_calib_df is not None else [],
            "sn_names": sn2_calib_df["SN"] if sn2_calib_df is not None else [],
            "magnitudes": sn2_calib_df["m0_Bi"] if sn2_calib_df is not None else [],
            "sigma": sn2_calib_df["sigma"] if sn2_calib_df is not None else [],
            # Alpha from Hubble flow fit
            "a_sn2": eq_data.get("a_sn2"),
            "a_sn2_err": eq_data.get("a_sn2_err"),
            "chisq_sn2_hf": eq_data.get("chisq_sn2_hf"),
            "ndof_sn2_hf": eq_data.get("ndof_sn2_hf"),
        } if sn2_calib_df is not None and len(sn2_calib_df) > 0 else {},

        # =====================================================================
        # TULLY-FISHER CALIBRATOR DETAILS
        # =====================================================================
        "tf_calib_details": {
            "n_tf": len(tf_calib_df) if tf_calib_df is not None else 0,
            "hosts": tf_calib_df["host"] if tf_calib_df is not None else [],
            "pgc": tf_calib_df["PGC"] if tf_calib_df is not None else [],
            "m": tf_calib_df["m"] if tf_calib_df is not None else [],
            "M": tf_calib_df["M"] if tf_calib_df is not None else [],
            "M_error": tf_calib_df["M_error"] if tf_calib_df is not None else [],
            # Alpha from Hubble flow fit
            "a_tf": eq_data.get("a_tf"),
            "a_tf_err": eq_data.get("a_tf_err"),
            "chisq_tf_hf": eq_data.get("chisq_tf_hf"),
            "ndof_tf_hf": eq_data.get("ndof_tf_hf"),
        } if tf_calib_df is not None and len(tf_calib_df) > 0 else {},

        # =====================================================================
        # SBF CALIBRATOR DETAILS
        # =====================================================================
        "sbf_calib_details": {
            "n_sbf": len(sbf_calib_df) if sbf_calib_df is not None else 0,
            "hosts": sbf_calib_df["host"] if sbf_calib_df is not None else [],
            "M110": sbf_calib_df["M110"] if sbf_calib_df is not None else [],
            "M110_error": sbf_calib_df["M110_error"] if sbf_calib_df is not None else [],
            "m110": sbf_calib_df["m110"] if sbf_calib_df is not None else [],
            "m110_error": sbf_calib_df["m110_error"] if sbf_calib_df is not None else [],
            # Alpha from Hubble flow fit
            "a_sbf": eq_data.get("a_sbf"),
            "a_sbf_err": eq_data.get("a_sbf_err"),
            "chisq_sbf_hf": eq_data.get("chisq_sbf_hf"),
            "ndof_sbf_hf": eq_data.get("ndof_sbf_hf"),
        } if sbf_calib_df is not None and len(sbf_calib_df) > 0 else {},

        # =====================================================================
        # COMA CLUSTER DETAILS
        # =====================================================================
        "coma_details": {
            "ncoma": eq_data.get("ncoma"),
            "ieq_coma_start": eq_data.get("ieq_coma_start"),
            "coma_sn_name": coma_df["coma_sn_name"] if coma_df is not None else [],
            "coma_sn_mag": coma_df["coma_sn_mag"] if coma_df is not None else [],
            "coma_sn_err": coma_df["coma_sn_err"] if coma_df is not None else [],
            "mu_coma_value": sol_results.get("mu_coma_value"),
            "mu_coma_error": sol_results.get("mu_coma_error")
        } if coma_df is not None else {},

        # =====================================================================
        # MEGAMASER DETAILS
        # =====================================================================
        "mm_details": {
            "nmm": eq_data.get("nmm"),
            "ieq_mm_start": eq_data.get("ieq_mm_start"),
            "mm_name": mm_df["name"] if mm_df is not None else [],
            "mm_logh_value": mm_df["logh_value"] if mm_df is not None else [],
            "mm_logh_error": mm_df["logh_error"] if mm_df is not None else [],
            "mm_v_total_error": mm_df["v_total_error"] if mm_df is not None else [],
            "vcorr_mm": read_dict.get("vcorr_mm"),
            "mm_v_corr": mm_df["v_corr"] if mm_df is not None else []
        } if mm_df is not None else {},

        # =====================================================================
        # EPM DETAILS
        # =====================================================================
        "epm_details": {
            "nepm": eq_data.get("nepm"),
            "ieq_epm_start": eq_data.get("ieq_epm_start"),
            "epm_sn": epm_df["SN"] if epm_df is not None else [],
            "epm_host": epm_df["Host"] if epm_df is not None else [],
            "epm_logh_value": epm_df["logh_value"] if epm_df is not None else [],
            "epm_logh_error": epm_df["logh_error"] if epm_df is not None else [],
            "vcorr_epm": read_dict.get("vcorr_epm")
        } if epm_df is not None else {},

        # =====================================================================
        # GROUPS DETAILS
        # =====================================================================
        "groups_details": {
            "n_groups": len(groups_df["group"].unique()) if groups_df is not None and len(groups_df) > 0 else 0,
            "ngroups_eq": eq_data.get("ngroups_eq"),
            "groups": groups_df["group"] if groups_df is not None else [],
            "group_hosts": groups_df["host"] if groups_df is not None else [],
            "group_sigma": groups_df["sigma"] if groups_df is not None else [],
            "group_param_index": eq_data.get("group_param_index")
        } if groups_df is not None and len(groups_df) > 0 else {},

        # =====================================================================
        # EQUALIZATION DROP DETAILS
        # =====================================================================
        "equalization_drop_details": {
            "removed_hosts_without_calibrators": removed_hosts_without_calibrators,
            "removed_sn1a_calib_without_hosts": removed_sn1a_calib_without_hosts,
            "removed_sn2_calib_without_hosts": removed_sn2_calib_without_hosts,
            "removed_tf_calib_without_hosts": removed_tf_calib_without_hosts,
            "removed_sbf_calib_without_hosts": removed_sbf_calib_without_hosts,
        }if removed_hosts_without_calibrators is not None or removed_sn1a_calib_without_hosts is not None or removed_sn2_calib_without_hosts is not None or removed_tf_calib_without_hosts is not None or removed_sbf_calib_without_hosts is not None else {},

        # =====================================================================
        # SOLUTION MATRICES
        # =====================================================================
        "params": sol_results.get("params"),
        "solmat": sol_results.get("solmat"),
        "invsolmat": sol_results.get("invsolmat"),
        "errmat": sol_results.get("errmat"),
        "covar": sol_results.get("covar"),
        "inv_covar": sol_results.get("inv_covar"),

        # =====================================================================
        # COVARIANCE DIAGNOSTICS
        # =====================================================================
        "covar_rank": sol_results.get("covar_rank"),
        "covar_cond": sol_results.get("covar_cond"),
        "covar_dim": sol_results.get("covar_dim"),
        "covar_rank_mc": eq_data.get("covar_rank_mc"),
        "ndof_full_mc": eq_data.get("ndof_full_mc"),
        "ndof_mc": eq_data.get("ndof_mc"),
        "covar_cond_mc": eq_data.get("covar_cond_mc"),
        "covar_dim_mc": eq_data.get("covar_dim_mc"),
    }

    return details




def make_json_serializable(obj):
    """Recursively convert numpy and pandas objects to JSON-serializable types.
    
    Handles nested dictionaries and lists, converting numpy arrays and pandas
    structures to native Python types.
    
    Args:
        obj (Any): Object to convert. Can be dict, list, np.ndarray, pd.Series,
            pd.DataFrame, numpy scalar, or native Python type.
    
    Returns:
        Any: JSON-serializable version of the input object.
            - np.ndarray and pd.Series -> list
            - pd.DataFrame -> dict with list values
            - numpy scalars -> Python scalars
            - dict/list -> recursively converted
            - pandas ExtensionArray -> list of Python scalars
            - other types -> unchanged
    """
    if isinstance(obj, dict):
        return {k: make_json_serializable(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [make_json_serializable(v) for v in obj]
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, pd.Series):
        return obj.tolist()
    elif isinstance(obj, pd.DataFrame):
        return obj.to_dict(orient="list")  # or orient="records" for list of row dicts
    elif isinstance(obj, (np.integer, np.floating)):
        return obj.item()
    elif hasattr(obj, '__array__'): # Handle pandas array types (StringArray, IntegerArray, etc.)
        return obj.tolist()
    else:
        return obj



def save_details_as_json(savefile, details):
    """Save detailed analysis results to JSON format.
    
    Converts the details dictionary to JSON-serializable format and writes to file.
    If the filename contains '.sav', it's replaced with '.json' extension.
    
    Args:
        savefile (str): Output filename. If contains '.sav', extension is changed to '.json'.
        details (dict): Details dictionary from generate_details().
    
    Example:
        >>> details = generate_details(...)
        >>> save_details_as_json('output.json', details)
    """
    # Convert details to format allowed by json
    details = make_json_serializable(details)
    # If it ends with ".sav" (or something like ".sav.Z"), replace everything after .sav with ".json"
    if ".sav" in savefile:
        base = savefile.split(".sav")[0]
        savefile = base + ".json"

    # Now save the dictionary
    with open(savefile, "w") as f:
        json.dump(details, f, indent=2)

    # Confirmation message
    config_reader.vprint(f"\nDetails saved to '{savefile}'")
