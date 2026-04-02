# H0DN: Distance Network Analysis Package

This package is a Python implementation developed for the H0DN collaboration and is based on the original IDL code written by Stefano Casertano. It constrains the Hubble constant (H0) in a distance network that combines various distance ladder approaches by taking into account the covariances of measurements. 

## Quick Start

### Installation

In the folder containing `pyproject.toml`, install in **editable mode**:

```bash
pip install -e .
```

### Basic Usage

Navigate to the `configs` folder and run the analysis:

```bash
cd configs
h0_constrainer config.ini
```

This performs a generalized least-squares fit using the configuration and data files in the data directory.

### Command-Line Options

```bash
# Use config.ini in current directory
h0_constrainer

# Use specific configuration file  
h0_constrainer my_config.ini

# Run multiple variants (systematic exploration)
h0_constrainer -v variants.ini
```

### Output

The tool prints results to console (when "verbose = False"):

```
H0=73.4988+/-0.8088 km/s/Mpc | Chi-squared=117.5597 | ndof=119 | Reduced Chi-squared=0.9879 | MB=-19.252+/-0.022
```

Optionally save detailed results and logs:

```ini
[DEFAULT]
logfile = output.txt      # Text summary
savefile = results.json    # Detailed results for analysis
```


## Configs

The `configs/` directory contains:

- **config.ini**: Working example configuration for the baseline analysis with extensive examples of alternative configs
- **variants.ini**: Multi-variant run (requires the base_variants_config.ini)
- **base_variants_config.ini**: A baseline configuration that the variants modify


## Data

The `data/` directory contains:
| Category | Files |
| --- | --- |
| Geometric distance measurements to anchors | `anchors.dat` |
| Host distances from primary distance indicators (cepheids, TRGB, etc.) | `host_data.dat`, `host_data_nomet.dat` |
| Calibrator magnitude data (SNe Ia, SBF, SN II, TF) | `sn1a_calib.dat`, `sn1a_calib_H.dat`, `sn1a_calib_J.dat`, `rung3_salt3.dat`, `sn1a_calib_bayessn.dat`, `sn1a_cal_csp_burns.dat`, `sbf_calib.dat`, `sn2_calib.dat`, `tf_calib.dat` |
| Hubble-flow data used to infer H0, and covariances | `hf_csp_lluis_v2.dat`, `hf_csp_lluis_v2_cov.dat`, `hf_salt3.dat`, `hf_salt3_cov.dat`, `sn1a_hf_H.dat`, `sn1a_hf_J.dat`, `sn1a_hf_bayessn.dat`, `sn1a_covar_bayessn.dat`, `sn1a_hf_pp.dat`, `sn1a_covar_pp.dat`, `sbf_hf.dat`, `sn2_hf.dat`, `tf_hf.dat` |
| Additional constraints and ties (megamasers, EPM, groups, coma) | `megamaser.dat`, `epm.dat`, `groups.dat`, `coma_sn.dat` |


Same files grouped by object/feataure type:
<table>
  <thead>
    <tr>
      <th>Object type</th>
      <th>Compilation</th>
      <th>Files</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>Anchors</td>
      <td>-</td>
      <td><code>anchors.dat</code></td>
    </tr>
    <tr>
      <td>Host distances</td>
      <td>-</td>
      <td><code>host_data.dat</code>, <code>host_data_nomet.dat</code></td>
    </tr>
    <tr>
      <td rowspan="7">SNe Ia</td>
      <td>PP</td>
      <td><code>sn1a_calib.dat</code>, <code>sn1a_hf_pp.dat</code>, <code>sn1a_covar_pp.dat</code></td>
    </tr>
    <tr>
      <td>BayesSN</td>
      <td><code>sn1a_calib_bayessn.dat</code>, <code>sn1a_hf_bayessn.dat</code>, <code>sn1a_covar_bayessn.dat</code></td>
    </tr>
    <tr>
      <td>CSP</td>
      <td><code>sn1a_cal_csp_burns.dat</code>, <code>hf_csp_lluis_v2.dat</code>, <code>hf_csp_lluis_v2_cov.dat</code></td>
    </tr>
    <tr>
      <td>SALT3</td>
      <td><code>hf_salt3.dat</code>, <code>hf_salt3_cov.dat</code>, <code>rung3_salt3.dat</code></td>
    </tr>
    <tr>
      <td>J</td>
      <td><code>sn1a_calib_J.dat</code>, <code>sn1a_hf_J.dat</code></td>
    </tr>
    <tr>
      <td>H</td>
      <td><code>sn1a_calib_H.dat</code>, <code>sn1a_hf_H.dat</code></td>
    </tr>
    <tr>
      <td>Coma</td>
      <td><code>coma_sn.dat</code></td>
    </tr>
    <tr>
      <td>SBF</td>
      <td>-</td>
      <td><code>sbf_calib.dat</code>, <code>sbf_hf.dat</code></td>
    </tr>
    <tr>
      <td>SN II</td>
      <td>-</td>
      <td><code>sn2_calib.dat</code>, <code>sn2_hf.dat</code>, <code>epm.dat</code> (not a cov file)</td>
    </tr>
    <tr>
      <td>Tully-Fisher</td>
      <td>-</td>
      <td><code>tf_calib.dat</code>, <code>tf_hf.dat</code></td>
    </tr>
    <tr>
      <td>Other</td>
      <td>-</td>
      <td><code>megamaser.dat</code>, <code>epm.dat</code>, <code>groups.dat</code>, <code>coma_sn.dat</code></td>
    </tr>
  </tbody>
</table>



## Package Structure

```
h0_constrainer/
├── __init__.py          # Package initialization
├── main.py              # Main workflow
├── cli.py               # Command-line interface
├── config_reader.py     # Configuration parsing
├── data_loader.py       # Data file loading
├── equations.py         # Equation system construction
├── solver.py            # Generalized least-squares solver
├── intercept.py         # Hubble flow intercept computation
├── logger.py            # Results organization
```


## Changelog

