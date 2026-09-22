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

from . import pipeline
def main(cf_file = None):
  pip = pipeline.H0DN_Pipeline()
  ini = pip.read_config(cf_file)
  data = pip.read_data(ini = ini)
  eq_data = pip.build_equations(ini = ini, data = data)
  sol_results = pip.solve(eq_data = eq_data)
  pip.log_results(ini = ini, sol_results = sol_results)
  pip.save_results(ini = ini, data = data, eq_data = eq_data, sol_results = sol_results)
  pip.print_results(ini = ini, sol_results = sol_results)
  
  print(pip.get_likelihood(eq_data, minimize_except=['M_B','H0']))
  class cosmo_class:
    def get_kz_function(self):
      return lambda z:1
    def get_cz_function(self):
      return lambda z:1
  cosmo = {"kz_function":(lambda z:1), "cz_function":(lambda z:1)}
  eq_data = pip.adjust_equations(cosmo)
  print(pip.get_likelihood(eq_data, minimize_except=['M_B','H0']))

if __name__ == "__main__":
    main()
