"""Data file loading and preprocessing module.

This module handles loading and preprocessing of all input data files for the
distance ladder analysis. It supports various file formats and distance indicators
including:
    - Geometric anchor distances [that calibrate the primary distance indicators (Cepheids, TRGB, etc.)]
    - Host galaxy distance measurements [inferred (mostly) from primary distance indicators]
    - Calibrator objects' magnitudes [calibrating secondary distance indicators (SNe Ia, SNe II, Tully-Fisher, SBF)]
    - Hubble flow samples for each secondary distance indicator
    - Additional constraints (megamasers, EPM SNe II, Coma cluster)
    - Galaxy groups for additional constraints on distances

Key Functions:
    load_anchors: Load geometric anchor calibrations.
    load_hosts: Load host galaxy distance measurements.
    load_calibrators: Load calibrator objects for various methods.
    load_hf: Load Hubble flow samples.
    load_groups: Load galaxy group membership for flow corrections.
    load_coma, load_mm, load_epm: Load additional distance constraints.
    drop_anchor_rows: Separate anchor calibrations from host measurements.
    compute_host_extras: Extract unique methods, anchors, and sources.
    uniquify: Create unique value arrays with index mappings.

Data Filtering:
    rung_equalize: Filter datasets to ensure consistent host coverage.
    should_exclude, should_include: Pattern-based inclusion/exclusion filtering.
    filter_groups: Remove groups with insufficient host or calibrator coverage.
"""

import os
import re
import io
import fnmatch
import warnings
import numpy as np
import pandas as pd
from scipy.constants import c as c_m_s
from h0_constrainer import config_reader

def load_rung_file(filepath, columns, skip_header_condition=None):
    """Load a whitespace-delimited data file with comment handling.
    
    Reads a file with '#' or ';' comments, optionally skipping header rows.
    
    Args:
        filepath (str): Path to the data file.
        columns (list): Column names for the DataFrame.
        skip_header_condition (callable, optional): Function that returns True
            for lines to skip (e.g., lambda line: line.startswith('#')).
    
    Returns:
        pd.DataFrame: Loaded data with specified column names.
    """
    with open(filepath, "r") as f:
        cleaned_data = [
            re.split(r"[#;]", line, maxsplit=1)[0].strip() 
            for line in f 
            if line.strip() and (skip_header_condition is None or not skip_header_condition(line))
        ]
    return pd.read_csv(io.StringIO("\n".join(cleaned_data)), sep=r"\s+", names=columns)




def load_hosts(datadir, filename):
    """Load host galaxy distance measurements.

    Loads the primary host data file containing distance modulus measurements
    from various methods (Cepheids, TRGB, etc.) anchored to geometric distances.
    Creates composite identifiers for method-anchor-source (MAS) and
    host-method-source (HMS) combinations.

    Args:
        datadir (str): Directory containing data files.
        filename (str): Host data filename. If empty/None, returns None.

    Returns:
        pd.DataFrame or None: DataFrame with columns:
            - host: Host galaxy name
            - mu_host: Distance modulus
            - mu_error: Distance modulus uncertainty
            - method: Distance indicator method (e.g., 'ceph_hst', 'trgb_jwst')
            - anchor: Geometric anchor (e.g., 'N4258', 'MW', 'LMC')
            - source: Data source/publication/collaboration
            - mas_name: Composite 'method&anchor&source' identifier
            - hms_name: Composite 'host&method&source' identifier
        Returns None if filename is empty.
    """
    if filename:
        path = os.path.join(datadir, filename)
        columns = ["host", "mu_host", "mu_error", "method", "anchor", "source","PGC"]
        df = load_rung_file(path, columns, skip_header_condition=lambda line: line.lower().startswith("host"))
        df["mas_name"] = df["method"] + "&" + df["anchor"] + "&" + df["source"]
        df["hms_name"] = df["host"] + "&" + df["method"] + "&" + df["source"]
        return df
    else:
        return None
    

def load_groups(datadir, filename):
    """Load galaxy group membership data.
    
    Groups are used to apply coherent flow corrections when multiple galaxies
    in the same physical group are observed; groups provide an additional 
    constraint on relative distances.
    
    Args:
        datadir (str): Directory containing data files.
        filename (str): Groups filename. If empty/None, returns None.
    
    Returns:
        pd.DataFrame or None: DataFrame with columns:
            - group: Group name
            - host: Host galaxy name (member of this group)
            - sigma: intrinsic distance modulus dispersion of the host with respect to the group
        Returns None if filename is empty.
    """
    if filename:
        path = os.path.join(datadir, filename)
        columns = ["group", "host", "sigma"]
        df = load_rung_file(
            path,
            columns,
            skip_header_condition=lambda line: line.strip().lower().startswith("groupname")
        )
        return df
    else:
        return None


def load_calibrators(datadir, filename, method):
    """Load calibrators for distance indicator methods.
    
    Loads calibrator objects/methods (SNe, TF galaxies, SBF observations) data.
    They should belong to hosts with measured distances (from primary distance indicators),
    this is handled by the filtering functions defined later.
    These calibrate the absolute magnitude parameter of the secondary distance indicator.
    
    Args:
        datadir (str): Directory containing data files.
        filename (str): Calibrator data filename. If empty/None, returns None.
        method (str): Method type - one of:
            - 'sn' or 'sn_IR': Type Ia supernovae
            - 'tf': Tully-Fisher relation
            - 'sbf': Surface Brightness Fluctuations
    
    Returns:
        pd.DataFrame or None: DataFrame with method-specific columns:
            For 'sn'/'sn_IR': host, SN, m0_Bi (magnitude), sigma
            For 'tf': host, PGC, m, M, M_error, mu_host, mu_error, anchor, source, method
            For 'sbf': host, M110, M110_error, m110, m110_error
        Returns None if filename is empty.
    
    Raises:
        ValueError: If method type is unknown.
    
    Note:
        IR mode ('sn_IR') handles optional J-band and H-band data formats.
    """
    if not filename:
        return None

    path = os.path.join(datadir, filename)

    if method == "sn":
        columns = ["host", "SN", "m0_Bi", "sigma"]
        df = load_rung_file(path, columns, skip_header_condition=lambda line: line.lower().startswith("host"))
        return df

    elif method == "tf":
        columns = ["host", "PGC", "m", "M", "M_error", "mu_host", "mu_error", "anchor", "source", "method"]
        df = load_rung_file(path, columns, skip_header_condition=lambda line: line.lower().startswith("host"))
        return df

    elif method == "sbf":
        columns = ["host", "M110", "M110_error", "m110", "m110_error"]
        df = load_rung_file(path, columns, skip_header_condition=lambda line: line.lower().startswith("host"))
        return df

    # --- IR calibrators: take first 4 *data* columns after optional index ---
    elif method == "sn_IR":
        # Read all, ignore '#' comments; whitespace-delimited.
        raw = pd.read_csv(path, sep=r"\s+", comment="#", header=None, engine="python")

        # Drop a possible bare header row (starts with 'gal', 'galaxy', or 'sn')
        first_tok = raw.iloc[:, 0].astype(str).str.strip().str.upper()
        raw = raw[~first_tok.isin(["GAL", "GALAXY", "SN"])]

        if raw.empty:
            return pd.DataFrame(columns=["host", "SN", "m0_Bi", "sigma"])

        # Detect optional leading index column (mostly numeric)
        idx0 = pd.to_numeric(raw.iloc[:, 0], errors="coerce")
        has_index = idx0.notna().mean() > 0.9
        off = 1 if has_index else 0

        need = off + 4
        if raw.shape[1] < need:
            raise ValueError(
                f"{filename}: found {raw.shape[1]} columns; need at least {need} "
                "(index? + host, SN, mag, emag)."
            )

        df = pd.DataFrame({
            "host":  raw.iloc[:, off + 0].astype(str),
            "SN":    raw.iloc[:, off + 1].astype(str),
            "m0_Bi": pd.to_numeric(raw.iloc[:, off + 2], errors="coerce"),
            "sigma": pd.to_numeric(raw.iloc[:, off + 3], errors="coerce"),
        })
        return df

    else:
        raise ValueError(f"Unknown method type: {method}")
    



def load_hf(datadir, filename, method, covar_filename=None):
    """
    Load "Hubble flow" (HF) rung data for various secondary distance indicators.

    This function reads Hubble flow (HF) data files for different methods
    (SNe Ia, SNe II, Tully-Fisher, SBF, SNe Ia IR), parses their columns
    according to the expected schema for each method, and returns a dataframe
    with standardized column names for downstream processing. Optionally,
    it can also load a covariance matrix for the measurements (if available).

    Args:
        datadir (str): Directory containing the input data files.
        filename (str): Data file name (relative to `datadir`).
        method (str): Indicator of the method/type of the Hubble flow rung.
            Accepted values:
                - "sn1a"    : Type Ia Supernovae
                - "sn2"     : Type II Supernovae
                - "tf"      : Tully-Fisher galaxies
                - "sbf"     : Surface Brightness Fluctuation
                - "sn1a_IR" : Type Ia Supernovae (IR bands, J or H)
        covar_filename (str, optional): Name of the covariance file (relative to `datadir`).
            If specified, a covariance matrix will be loaded and returned alongside the data.

    Returns:
        (df, np.ndarray):
            - If `covar_filename` is provided, returns a tuple of dataframe and the covariance matrix.
            - If `covar_filename` is not provided, returns a tuple of dataframe and None.

    Raises:
        ValueError: If an unknown method is specified or if the file format does not match expectations.

    Notes:
        - The columns read from the files and written to the dataframe are method-dependent.
        - For "sn1a_IR", if more than one magnitude band is present (10 columns), both are loaded ("mb2", "mb_err2").
        - For SBF/Tully-Fisher/SN1a files, comment lines with '#' at start are skipped.
        - For covariance files, the first line denoting size is skipped, covariance returned as a square matrix.

    Example:
        >>> df, cov = load_hf('data', 'sn1a_hf.dat', 'sn1a')
        >>> df, cov = load_hf('data', 'sn1a_hf.dat', 'sn1a', 'sn1a_cov.dat')
    """
    if not filename:
        return None, None

    path = os.path.join(datadir, filename)

    # --- Standard methods ---
    if method == "sn1a":
        columns = ["name", "mb", "mb_err", "zhel", "zcmb", "vp", "vp_2mpp_sdss_6df", "vp_2mrs", "vp_2mpp"]
        df = load_rung_file(path, columns, skip_header_condition=lambda line: line.strip().startswith("#"))

    elif method == "sn2":
        columns = ["name", "z_corr", "m0_i", "em0_i"]
        df = load_rung_file(path, columns, skip_header_condition=lambda line: line.strip().startswith("#"))

    elif method == "tf":
        columns = ["name", "PCG", "mu", "mu_err", "M", "m", "sigma", "v_cmb", "v_2mpp", "v_2mrs", "v_2mpp_6df_sdss"]
        df = load_rung_file(path, columns, skip_header_condition=lambda line: line.strip().startswith("#"))

    elif method == "sbf":  # SBF
        columns = ["Name", "M110", "e(M110)", "m110", "e(m110)", "vgrp", "v2m++"]
        df = load_rung_file(path, columns, skip_header_condition=lambda line: line.strip().startswith("#"))

    # --- For SN Ia IR HF files (J/H order agnostic) ---
    elif method == "sn1a_IR":
        # Read everything (whitespace-delimited), ignore lines starting with '#'.
        raw = pd.read_csv(path, sep=r"\s+", comment="#", header=None, engine="python")

        # Drop a possible bare header row (e.g., "SN GALAXY ...")
        mask_header = raw.iloc[:, 0].astype(str).str.strip().str.upper().eq("SN")
        raw = raw[~mask_header]

        # We expect at least 8 columns:
        # 0: name, 1: galaxy, 2-3: first mag pair (J or H), 4: zhel, 5: zcmb,
        # 6: ezcmb, 7: zcorr, 8-9: second mag pair (optional)
        if raw.shape[1] < 8:
            raise ValueError(
                f"sn1a_IR file has {raw.shape[1]} columns; expected at least 8 "
                "(name, galaxy, mag1, emag1, zhel, zcmb, ezcmb, zcorr, [mag2, emag2])."
            )

        # Build the output with neutral names; keep both mag pairs if present.
        df = pd.DataFrame({
            "name":   raw.iloc[:, 0].astype(str),
            "galaxy": raw.iloc[:, 1].astype(str),
            "mb":   pd.to_numeric(raw.iloc[:, 2], errors="coerce"),
            "mb_err":  pd.to_numeric(raw.iloc[:, 3], errors="coerce"),
            "zhel":   pd.to_numeric(raw.iloc[:, 4], errors="coerce"),
            "zcmb":   pd.to_numeric(raw.iloc[:, 5], errors="coerce"),
            "ezcmb":  pd.to_numeric(raw.iloc[:, 6], errors="coerce"),
            "zcorr":  pd.to_numeric(raw.iloc[:, 7], errors="coerce"),
        })
        if raw.shape[1] >= 10:
            df["mb2"]  = pd.to_numeric(raw.iloc[:, 8], errors="coerce")
            df["mb_err2"] = pd.to_numeric(raw.iloc[:, 9], errors="coerce")
        # print(df)

    else:
        raise ValueError(f"Unknown method type: {method}")

    # --- Optional covariance handling ---
    if covar_filename:
        cov_path = os.path.join(datadir, covar_filename)
        cov = np.loadtxt(cov_path)[1:]  # skip the first line (it redundantly denotes n for nxn cov matrix)
        lcov = int(np.sqrt(len(cov)))
        return df, cov.reshape(lcov, lcov)
    else:
        return df, None


def drop_anchor_rows(host_df):
    """Filters the host data file to remove anchor rows (host=anchor) and create MAS reference arrays.

    This function identifies anchor rows in the host data,
    extracts their reference values, removes them from the working dataset,
    and creates Method-Anchor-Source (MAS) tracking arrays to later and index mappings to
    update the covariance matrix later in equations.py.

    Args:
        host_df (pd.DataFrame): Host data including anchor rows.

    Returns:
        tuple: (host_df_filtered, mas_ref_error, mas_ref_value, mas_ref_anchor,
                mas_ref_method, mas_ref_source, n_mas, mas, mas_index, mas_error)
            - host_df_filtered: Host data with anchor rows removed.
            - mas_ref_error: Uncertainties for each anchor MAS combination.
            - mas_ref_value: Distance moduli for each anchor MAS combination.
            - mas_ref_anchor: Anchor names for each reference.
            - mas_ref_method: Method names for each reference.
            - mas_ref_source: Source names for each reference.
            - n_mas: Number of unique MAS combinations.
            - mas: Array of unique MAS names.
            - mas_index: Index mapping each host row to its MAS group.
            - mas_error: Uncertainty for each unique MAS (matched to anchors).

    Note:
        Anchor list: {'N4258', 'MW', 'LMC', 'SMC', 'M31', 'ALL'}
        MAS uncertainties are matched by comparing MAS names constructed as
        'method&anchor&source' between host rows and anchor references.
    """
    anchor_list = {"N4258", "MW", "LMC", "SMC", "M31", "ALL"}
    is_anchor = host_df["host"].isin(anchor_list)

    # Save MAS reference values from anchor rows
    mas_ref_error = host_df.loc[is_anchor, "mu_error"].values
    mas_ref_value = host_df.loc[is_anchor, "mu_host"].values
    mas_ref_anchor = host_df.loc[is_anchor, "anchor"].values
    mas_ref_method = host_df.loc[is_anchor, "method"].values
    mas_ref_source = host_df.loc[is_anchor, "source"].values

    # Drop the anchor rows and reset the index
    host_df_filtered = host_df.loc[~is_anchor].reset_index(drop=True)

    # Extract unique MAS groups from host data
    mas, mas_index = uniquify(host_df_filtered["mas_name"].values)
    n_mas = len(mas)

    # Construct MAS reference names from the anchor rows that were dropped
    if mas_ref_method.size > 0 and mas_ref_anchor.size > 0 and mas_ref_source.size > 0:
        mas_ref_name = np.array(mas_ref_method) + "&" + np.array(mas_ref_anchor) + "&" + np.array(mas_ref_source)
    else:
        mas_ref_name = np.array([])

    # Initialize MAS uncertainties with a default value
    mas_error = np.ones(n_mas, dtype=float)

    # For each unique MAS group, update the uncertainty if a matching reference exists
    for i, mas_val in enumerate(mas):
        if mas_ref_name.size > 0:
            match_idx = np.where(mas_ref_name == mas_val)[0]
            if len(match_idx) > 0:
                if match_idx[0] < len(mas_ref_error):
                    mas_error[i] = mas_ref_error[match_idx[0]]
                else:
                    print(f"Error: match_idx {match_idx[0]} out of bounds for mas_ref_error")
            else:
                config_reader.wprint(f"=== WARNING === No reference error found for {mas_val}. Using default uncertainty (1.0).")
        else:
            config_reader.vprint(f"No MAS reference data available for '{mas_val}'; using default uncertainties.")

    return (host_df_filtered,
            mas_ref_error, mas_ref_value, mas_ref_anchor,
            mas_ref_method, mas_ref_source, n_mas, mas, mas_index, mas_error)



def load_anchors(datadir, anchorfile="anchors.dat", default_anchors=["N4258", "MW", "LMC", "SMC"], default_mu_sigma=[0.032, 0.001, 0.024, 0.032]):
    """Load geometric anchor distance data.
    
    Reads geometric anchor distance moduli and their uncertainties. These
    anchors provide the absolute distance scale for the distance ladder.
    Falls back to default values if file is missing or malformed.
    
    Args:
        datadir (str): Directory containing data files.
        anchorfile (str, optional): Anchor data filename. Defaults to 'anchors.dat'.
        default_anchors (list, optional): Default anchor names.
        default_mu_sigma (list, optional): Default geometric uncertainties.
    
    Returns:
        tuple: (anchor_data, anchor_names, anchor_g_mu_sigma)
            - anchor_data (pd.DataFrame or None): Full anchor data table.
            - anchor_names (np.ndarray): Array of anchor names.
            - anchor_g_mu_sigma (np.ndarray): Array of geometric uncertainties.
    
    Note:
        Expected file format: whitespace-delimited with columns:
        anchor, mu_value, mu_sigma
        Lines starting with ';' are treated as comments.
    """
    anchors_file = os.path.join(datadir, anchorfile)
    if os.path.exists(anchors_file):
        cleaned_lines = []
        with open(anchors_file, "r") as file:
            for line in file:
                clean_line = re.split(r";", line, maxsplit=1)[0].strip()
                if clean_line:
                    cleaned_lines.append(clean_line)
        try:
            anchor_data = pd.read_csv(io.StringIO("\n".join(cleaned_lines)), sep=r"\s+", 
                                      names=["anchor", "mu_value", "mu_sigma"], skiprows=1)
            anchor_data["anchor"] = anchor_data["anchor"].astype(str).str.strip()
            anchor_names = anchor_data["anchor"].values
            anchor_g_mu_sigma = anchor_data["mu_sigma"].values.astype(float)
            return anchor_data, anchor_names, anchor_g_mu_sigma
        except Exception as e:
            print(f"Error loading anchors.dat: {e}")
            return None, np.array(default_anchors), np.array(default_mu_sigma)
    else:
        print(f"Error loading anchors.dat: {anchors_file}, defaulting to 'anchors.dat'")
        return None, np.array(default_anchors), np.array(default_mu_sigma)

def compute_host_extras(host_df, anchor_names, anchor_g_mu_sigma):
    """Extract unique methods, anchors, sources and assign uncertainties.

    Analyzes the host data to identify all unique methods (cepheids_hst, trgb_jwst etc.), 
    geometric anchors (e.g., NGC 4258, LMC), and data sources (R22, CCHP etc.). 
    Creates index mappings and assigns appropriate uncertainties for
    anchors and host-method-source combinations.

    Args:
        host_df (pd.DataFrame): Host distance data.
        anchor_names (np.ndarray): Array of anchor names from load_anchors().
        anchor_g_mu_sigma (np.ndarray): Array of anchor uncertainties.

    Returns:
        tuple: (methods, method_index, anchors, anchor_index, sources, source_index,
                hms, hms_index, mu_hms_error, mu_anchor_error)
            - methods: Unique method names.
            - method_index: Index mapping each host row to its method.
            - anchors: Unique anchor names.
            - anchor_index: Index mapping each host row to its anchor.
            - sources: Unique source names.
            - source_index: Index mapping each host row to its source.
            - hms: Unique host-method-source combinations.
            - hms_index: Index mapping each host row to its HMS group.
            - mu_hms_error: Uncertainty for each HMS combination.
            - mu_anchor_error: Geometric uncertainty for each anchor.

    Warnings:
        Issues warnings if HMS combinations have inconsistent reported errors.
    """
    methods, method_index = uniquify(host_df["method"].values)
    anchors, anchor_index = uniquify(host_df["anchor"].values)
    sources, source_index = uniquify(host_df["source"].values)

    # Assign uncertainties to anchors based on geometric distance uncertainties
    mu_anchor_error = np.ones(len(anchors))
    for i, anchor in enumerate(anchors):
        anchor_cleaned = str(anchor).strip()
        # Find matching anchor index
        match_idx = np.where(np.array([s.strip() for s in anchor_names.astype(str)]) == anchor_cleaned)[0]
        if len(match_idx) > 0:
            mu_anchor_error[i] = anchor_g_mu_sigma[match_idx[0]]

    # Assign uncertainties for host-method-source combinations
    hms, hms_index = uniquify(host_df["hms_name"].values)
    mu_hms_error = np.zeros(len(hms))
    for i, name in enumerate(hms):
        wh = np.where(hms_index == i)[0]
        hms_errors = host_df["mu_error"].values[wh]
        if np.allclose(hms_errors, hms_errors[0]):
            mu_hms_error[i] = hms_errors[0]
        else:
            mu_hms_error[i] = np.min(hms_errors)
            config_reader.wprint(f"\n=== WARNING === Inconsistent mu_error values for HMS '{name}', using minimum = {mu_hms_error[i]:.4f}")

    return methods, method_index, anchors, anchor_index, sources, source_index, hms, hms_index, mu_hms_error, mu_anchor_error

def load_coma(coma_file):
    """
    Load SN/Coma cluster constraint data.

    Reads a file of SN/Coma measurements, returning a DataFrame with columns:
        - coma_sn_name: SN name or ID
        - coma_sn_mag:  SN apparent magnitude
        - coma_sn_err:  Magnitude uncertainty

    Returns None if the file is empty or not specified.
    """
    if coma_file and coma_file.lower() != "none" and coma_file.strip():
        with open(coma_file, "r") as file:
            cleaned_lines = [re.split(r"#", line, maxsplit=1)[0].strip() 
                             for line in file if line.strip()]
        coma_columns = ["coma_sn_name", "coma_sn_mag", "coma_sn_err"]
        df = pd.read_csv(io.StringIO("\n".join(cleaned_lines)), sep=r"\s+", names=coma_columns)
        return df
    else:
        return None

def load_mm(mm_file, q0, j0, vpec_error, vcorr_mm="2M++", mm_dz=True):
    """
    Loads the mm constraints file and processes it.
    
    If mm_dz is True (default value), then the file is expected to contain 
    distance and redshift information. In that case, additional columns are computed:
      - distance (used to compute the dist modulus)
      - distance_error (from d_plus and d_minus)
      - v_corr (selecting the velocity correction column based on vcorr_mm)
      - z (redshift, computed as v_corr / c)
      - cz_term, v_total_error, cz_error (ancillary quantities)
      - logh_value, logh_error (final computed quantities)
      
    If mm_dz is False, then it is assumed that H0 and its error are provided directly.
    """
    if mm_file and mm_file.lower() != "none" and mm_file.strip():
        with open(mm_file, "r") as file:
            cleaned_lines = [
                re.split(r"#", line, maxsplit=1)[0].strip() 
                for line in file if line.strip()
            ]
        # For the distance-redshift case, we expect the following columns:
        column_names = [
            "name", "distance", "d_minus", "d_plus", "v_obs", "v_error",
            "v_corr_group", "v_corr_2M++", "v_corr_CF3", "v_corr_M2000",
            "pec_group", "pec_2M++", "pec_CF3", "pec_M2000"
        ]
        df = pd.read_csv(io.StringIO("\n".join(cleaned_lines)), sep=r"\s+", names=column_names)
        
        if mm_dz:
            # Compute the distance error from d_plus and d_minus
            df["distance_error"] = 0.5 * (df["d_plus"].abs() + df["d_minus"].abs())
            # Compute the distance modulus and its error
            df["mu"] = 5 * np.log10(df["distance"]) + 25
            df["mu_error"] = 5 * (df["distance_error"] / df["distance"]) / np.log(10)
            # Select the velocity correction column based on vcorr_mm
            vcorr_map = {
                "GROUP": "v_corr_group",
                "2M++": "v_corr_2M++",
                "CF3": "v_corr_CF3",
                "M2000": "v_corr_M2000"
            }
            vcorr_column = vcorr_map.get(vcorr_mm.upper(), "v_corr_2M++")
            if vcorr_column in df.columns:
                df["v_corr"] = df[vcorr_column]
            else:
                config_reader.wprint(f"=== WARNING === {vcorr_column} not found in mm file. Using v_corr_2M++ as default.")
                df["v_corr"] = df["v_corr_2M++"]
            # Compute redshift (v_corr/c)
            df["z"] = df["v_corr"] / (c_m_s/1000.0)
            # Compute the velocity correction term
            df["cz_term"] = np.log10(df["v_corr"] * (1.0 - 0.5*(q0+3) * df["z"] + 1./6.*(11.0 + 7.0 * q0 + 3.0 * q0**2 - j0) * df["z"]**2))
            # Compute total velocity error (combining measurement and peculiar velocity errors)
            df["v_total_error"] = np.sqrt(df["v_error"]**2 + vpec_error**2)
            df["cz_error"] = df["v_total_error"] / df["v_corr"] / np.log(10)
            # Compute log10(H0) value and its error
            df["logh_value"] = df["cz_term"] - 0.2 * df["mu"] + 5
            df["logh_error"] = np.sqrt((0.2 * df["mu_error"])**2 + df["cz_error"]**2)
        else:
            # Direct case: H0 and its error are provided directly.
            df["logh_value"] = np.log10(df["h0_value"])
            df["logh_error"] = df["h0_error"] / (df["h0_value"] * np.log(10))
        return df
    else:
        return None
    


def load_epm(epm_file, q0, j0, vpec_error, vcorr_epm="2M++", epm_dz=True):
    """
    Loads the epm constraints file and processes it.
    
    If epm_dz is True (default value), then the file is expected to contain 
    distance and redshift information. In that case, additional columns are computed:
      - distance, distance_error (used to compute the dist modulus)
      - v_corr (selecting the velocity correction column based on vcorr_epm)
      - z (redshift, computed as v_corr / c)
      - cz_term, v_total_error, cz_error (ancillary quantities)
      - logh_value, logh_error (final computed quantities)
      
    If epm_dz is False, then it is assumed that H0 and its error are provided directly.
    """
    if epm_file and epm_file.lower() != "none" and epm_file.strip():
        with open(epm_file, "r") as file:
            cleaned_lines = [
                re.split(r"#", line, maxsplit=1)[0].strip()
                for line in file if line.strip()
            ]
        # For the distance-redshift case, we expect the following columns:
        column_names = [
            "No", "Host", "SN", "distance", "distance_error",
            "z_cmb", "v_pec_2M++", "z_cosmo", "Source"
        ]
        df = pd.read_csv(io.StringIO("\n".join(cleaned_lines)), sep=r"\s+", names=column_names)

        if epm_dz:
            # Compute the distance modulus and its error
            df["mu"] = 5 * np.log10(df["distance"]) + 25
            df["mu_error"] = 5 * (df["distance_error"] / df["distance"]) / np.log(10)
            # Select the redshift column based on vcorr_epm
            vcorr_map = {
                "CMB": "z_cmb",
                "2M++": "z_cosmo"
            }
            z_column = vcorr_map.get(vcorr_epm.upper(), "z_cosmo")
            if z_column in df.columns:
                df["v_corr"] = df[z_column] * (c_m_s / 1000.0)
            else:
                config_reader.wprint(f"=== WARNING === {z_column} not found in epm file. Using z_cosmo as default.")
                df["v_corr"] = df["z_cosmo"] * (c_m_s / 1000.0)
            # Compute redshift (v_corr / c)
            df["z"] = df["v_corr"] / (c_m_s / 1000.0)
            # Compute the velocity correction term
            df["cz_term"] = np.log10(
                df["v_corr"] * (1.0 + (1.0 - q0) * df["z"] / 2.0 - (1.0 - q0 - 3.0 * q0**2.0 + j0)/6.0 * df["z"]**2.0)
            )
            # Compute total velocity error (peculiar + nominal 1 km/s measurement)
            df["v_error"] = 1.0
            df["v_total_error"] = np.sqrt(df["v_error"]**2 + vpec_error**2)
            df["cz_error"] = df["v_total_error"] / df["v_corr"] / np.log(10)
            # Compute log10(H0) value and its error
            df["logh_value"] = df["cz_term"] - 0.2 * df["mu"] + 5
            df["logh_error"] = np.sqrt((0.2 * df["mu_error"])**2 + df["cz_error"]**2)
        else:
            # Direct case: H0 and its error are provided directly.
            df["logh_value"] = np.log10(df["h0_value"])
            df["logh_error"] = df["h0_error"] / (df["h0_value"] * np.log(10))
        return df
    else:
        return None




def uniquify(array):
    """Create unique values array with index mapping.
    
    Takes an array of (possibly repeated) values and returns both the unique
    values and an index array that maps each original element to its unique value.
    
    Args:
        array (array-like): Input array with potentially repeated values.
    
    Returns:
        tuple: (unique_vals, index_map)
            - unique_vals (np.ndarray): Sorted unique values from input.
            - index_map (np.ndarray): For each element in original array,
              the index into unique_vals. Same shape as input.
    
    Example:
        >>> vals, idx = uniquify(['A', 'B', 'A', 'C', 'B'])
        >>> vals
        array(['A', 'B', 'C'])
        >>> idx
        array([0, 1, 0, 2, 1])
    """
    array = np.array(array)
    sorted_indices = np.argsort(array)
    sorted_array = array[sorted_indices]
    unique_vals, unique_pos = np.unique(sorted_array, return_index=True)
    index_map = np.zeros_like(sorted_indices, dtype=int)
    index_map[sorted_indices] = np.searchsorted(unique_vals, sorted_array)
    return unique_vals, index_map




# Could be useful for future data formats but not used yet:
# def matches_condition(value, pattern):
#     """Check if the value matches the wildcard pattern (case-sensitive)."""
#     return fnmatch.fnmatchcase(str(value), pattern)

def matches_condition(value, pattern):
    """Check if value matches a wildcard pattern (case-insensitive).
    
    Args:
        value: Value to test.
        pattern (str): Wildcard pattern (* and ? supported).
    
    Returns:
        bool: True if value matches pattern.
    """
    return fnmatch.fnmatch(str(value).lower(), pattern.lower())


def should_exclude(row, exclude_list, columns):
    """Check if a row should be excluded based on pattern list.
    
    Args:
        row (pd.Series): DataFrame row to test.
        exclude_list (list): List of exclusion patterns. Multiple conditions
            within a pattern are joined with '&'. Column-specific patterns use
            'column:pattern' syntax. See example below.
        columns (list): Column names to search if no column specified.
    
    Returns:
        bool: True if row matches any exclusion pattern.
    
    Example:
        >>> should_exclude(row, ['N4258&ceph*', 'source:R22', 'host:N4258&method:ceph_hst'], ['host', 'method'])
    """
    for exclude_entry in exclude_list:
        conditions = exclude_entry.split("&")
        if all(
            (
                matches_condition(row[term.split(":")[0]], term.split(":")[1])
                if ":" in term and term.split(":")[0] in row
                else any(matches_condition(row[col], term) for col in columns)
            )
            for term in conditions
        ):
            return True
    return False

def should_include(row, inclusion_list, columns):
    """Check if a row should be included based on pattern list.
    
    Args:
        row (pd.Series): DataFrame row to test.
        inclusion_list (list): List of inclusion patterns (same syntax as should_exclude).
        columns (list): Column names to search if no column specified.
    
    Returns:
        bool: True if row matches any inclusion pattern.
    """
    for include_entry in inclusion_list:
        conditions = include_entry.split("&")
        if all(
            (
                matches_condition(row[term.split(":")[0]], term.split(":")[1])
                if ":" in term and term.split(":")[0] in row
                else any(matches_condition(row[col], term) for col in columns)
            )
            for term in conditions
        ):
            return True
    return False



def rung_equalize(df_ref, df_target, reverse=False, df_groups=None):
    """Filter datasets to ensure consistent host coverage across distance ladder rungs.
    
    Ensures that calibrator objects only exist for hosts that have distance measurements,
    and vice versa. In forward mode, removes rows from the target with hosts that are not the reference DataFrame. In reverse
    mode, removes rows from the reference with hosts that are not the target DataFrame (keeping anchor rows).

    Args:
        df_ref (pd.DataFrame): Reference DataFrame providing allowed host names (typically host_df).
        df_target (pd.DataFrame): DataFrame to filter (typically calibrator DataFrame).
        reverse (bool, optional): If True, filter df_ref based on df_target instead
            of filtering df_target based on df_ref. Defaults to False.
        df_groups (pd.DataFrame, optional): Groups DataFrame. If provided, group hosts
            are also treated as valid. Defaults to None.

    Returns:
        tuple:
            - pd.DataFrame: Filtered DataFrame with only hosts that have coverage in both
                datasets (plus anchors in reverse mode).
            - list or None: List of removed entries that were dropped during filtering,
                or None if none were removed.

    Note:
        Forward mode: Used to remove calibrators for hosts without distance measurements.
        Reverse mode: Used to remove host measurements without calibrators (keeping anchors).
    """
    anchor_list = {"N4258", "MW", "LMC", "SMC", "M31", "ALL"}

    # Collect extra hosts from groups if provided
    extra_hosts = set()
    if df_groups is not None and len(df_groups) > 0:
        extra_hosts = set(df_groups["host"].astype(str))

    if not reverse:
        valid_hosts = set(df_ref["host"]) | extra_hosts
        removed_hosts = sorted(set(df_target["host"]) - valid_hosts)
        filtered_df = df_target[df_target["host"].isin(valid_hosts)]
        config_reader.vprint(
            f" calibrator hosts: {len(filtered_df)} rows after removing objects with unknown host distances. "
            f"Removed hosts: {', '.join(removed_hosts) if removed_hosts else 'None'}"
        ) #This vprint makes sense when taken with the previous print statement in main.py
    else:
        valid_hosts = set(df_target["host"]) | extra_hosts
        is_valid = df_ref["host"].isin(valid_hosts)
        is_anchor = df_ref["host"].isin(anchor_list)
        keep_mask = is_valid | is_anchor
        removed_hosts = sorted(set(df_ref["host"][~keep_mask]))
        filtered_df = df_ref[keep_mask]
        config_reader.vprint(
            f"Host distance measurements file: {len(filtered_df)} rows (including anchors) "
            f"after removing {len(removed_hosts)} hosts "
            f"with no calibrator objects."
        )

    return filtered_df, removed_hosts if removed_hosts else None


def filter_groups(host_df, groups_df,
                  sn1a_calib_df=None,
                  sn2_calib_df=None,
                  tf_calib_df=None,
                  sbf_calib_df=None):
    """Filter groups based on host and calibrator coverage.

    Applies three filtering rules:
    1) Drop a group if NONE of its hosts are present in host_df['host'].
    2) Drop a group if NONE of its hosts are present in ANY calibrator DataFrame['host'].
    3) Drop host-rows from groups_df if that host appears multiple times in host_df['host'].

    Args:
        host_df (pd.DataFrame): Host distance measurements.
        groups_df (pd.DataFrame): Galaxy group membership data.
        sn1a_calib_df (pd.DataFrame, optional): SNe Ia calibrators.
        sn2_calib_df (pd.DataFrame, optional): SNe II calibrators.
        tf_calib_df (pd.DataFrame, optional): Tully-Fisher calibrators.
        sbf_calib_df (pd.DataFrame, optional): SBF calibrators.

    Returns:
        pd.DataFrame: Filtered groups DataFrame.
    """
    groups_df["group"] = groups_df["group"].astype(str)
    groups_df["host"] = groups_df["host"].astype(str)

    available_hosts = set(host_df["host"].astype(str))

    # Collect calibrator hosts from all provided DataFrames
    calib_hosts = set()
    for df in [sn1a_calib_df, sn2_calib_df, tf_calib_df, sbf_calib_df]:
        if df is not None and len(df) > 0:
            calib_hosts.update(df["host"].astype(str).tolist())

    # Rule 1: Drop groups with no overlap to host_df hosts
    grp_to_hosts = groups_df.groupby("group")["host"].apply(set)
    to_drop_no_hosts = [grp for grp, hosts in grp_to_hosts.items() if hosts.isdisjoint(available_hosts)]
    for grp in sorted(to_drop_no_hosts):
        config_reader.vprint(f"Dropping group {grp}: no host belonging to this group exists in host_data file.")
    groups_df = groups_df[~groups_df["group"].isin(to_drop_no_hosts)].reset_index(drop=True)

    if len(groups_df) == 0:
        return groups_df

    # Rule 2: Drop groups with no overlap to calibrator hosts
    if calib_hosts:
        grp_to_hosts = groups_df.groupby("group")["host"].apply(set)
        to_drop_no_calib = [grp for grp, hosts in grp_to_hosts.items() if hosts.isdisjoint(calib_hosts)]
        for grp in sorted(to_drop_no_calib):
            config_reader.vprint(f"Dropping group {grp}: no host belonging to this group exists in any calib file.")
        groups_df = groups_df[~groups_df["group"].isin(to_drop_no_calib)].reset_index(drop=True)

    if len(groups_df) == 0:
        return groups_df

    # Rule 3: Drop host-rows where host has multiple entries in host_df
    counts = host_df["host"].astype(str).value_counts()
    dup_hosts = set(counts[counts > 1].index.tolist())

    if dup_hosts:
        mask_dup = groups_df["host"].isin(dup_hosts)
        if mask_dup.any():
            for grp, host in groups_df.loc[mask_dup, ["group", "host"]].itertuples(index=False):
                config_reader.vprint(f"Warning: host {host} not considered as a member of group {grp}; "
                                     f"multiple entries of this host in the host_data file.")
            groups_df = groups_df.loc[~mask_dup].reset_index(drop=True)

    return groups_df


