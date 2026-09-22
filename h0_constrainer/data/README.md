Folder with data files for both IDL and Python versions. See instructions in the idlcode folder for information on how to direct the code to look for data files here.

Several of the data files here are used in variants.  Please check the code in idlcode/vars.pro or the instructions in h0_constrainer/config to identify which files are used in which variants.


### Summary of files in this folder
Data sets used in the paper's baseline analysis are marked in bold below. Note: in the baseline configuration, `host_data.dat` excludes rows that contain any of 'jagb', 'mira', 'b25', or 'SMC'.
| Category | Files |
| --- | --- |
| Geometric distance measurements to anchors | `anchors.dat` |
| Host distances from primary distance indicators (cepheids, TRGB, etc.) | **`host_data.dat`**, `host_data_nomet.dat` |
| Calibrator magnitude data (SNe Ia, SBF, SN II, TF) | **`sn1a_calib.dat`**, `sn1a_calib_H.dat`, `sn1a_calib_J.dat`, `rung3_salt3.dat`, `sn1a_calib_bayessn.dat`, `sn1a_cal_csp_burns.dat`, **`sbf_calib.dat`**, `sn2_calib.dat`, `tf_calib.dat` |
| Hubble-flow data used to infer H0, and covariances | `hf_csp_lluis_v2.dat`, `hf_csp_lluis_v2_cov.dat`, `hf_salt3.dat`, `hf_salt3_cov.dat`, `sn1a_hf_H.dat`, `sn1a_hf_J.dat`, `sn1a_hf_bayessn.dat`, `sn1a_covar_bayessn.dat`, **`sn1a_hf_pp.dat`**, **`sn1a_covar_pp.dat`**, **`sbf_hf.dat`**, `sn2_hf.dat`, `tf_hf.dat` |
| Additional constraints and ties (megamasers, EPM, groups, coma) | **`megamaser.dat`**, `epm.dat`, **`groups.dat`**, `coma_sn.dat` |


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
      <td><strong><code>host_data.dat</code></strong>, <code>host_data_nomet.dat</code></td>
    </tr>
    <tr>
      <td rowspan="7">SNe Ia</td>
      <td>PP</td>
      <td><strong><code>sn1a_calib.dat</code></strong>, <strong><code>sn1a_hf_pp.dat</code></strong>, <strong><code>sn1a_covar_pp.dat</code></strong></td>
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
      <td><strong><code>sbf_calib.dat</code></strong>, <strong><code>sbf_hf.dat</code></strong></td>
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
      <td><strong><code>megamaser.dat</code></strong>, <code>epm.dat</code>, <strong><code>groups.dat</code></strong>, <code>coma_sn.dat</code></td>
    </tr>
  </tbody>
</table>

