from cobaya.likelihood import Likelihood
from h0_constrainer import pipeline
from scipy.interpolate import CubicSpline
import numpy as np
import os

from pathlib import Path
import h0_constrainer

_DEFAULT_CONFIG_DIR = str(Path(h0_constrainer.__file__).resolve().parents[0] / "configs" / "base_config.ini")
class H0DN_like(Likelihood):
    zinterp : float = np.linspace(1e-5, 3, num=1000)
    cf_file : str = _DEFAULT_CONFIG_DIR
    minimize_except : list[str] = ["M_B","H0"]
    icov  = None
    #modify_config : dict[str] = {"sn1a_hf_file" : "", "sn1a_hf_cov_file":"", "alpha_sn1a":0.714, "alpha_sn1a_error":1000.0}
    modify_config : dict[str] = {}

    def initialize(self):
        """Called once during sampler initialization."""
        self.pip = pipeline.H0DN_Pipeline()
        self.pip.read_config(os.path.abspath(self.cf_file))
        self.pip.ini.update(self.modify_config)
        self.pip.read_data()
        self.pip.build_equations()
        self.pip.solve() # To get inverse covariance matrix, to make subsequent likelihood constructions faster
        # Get likelihood at example point to get useful info like reasonable bestift, covmat, ...
        self.icov = self.pip.sol_results['inv_covar']
        pars, lkl, bf, cov = self.pip.get_likelihood(inv_covar = self.icov, minimize_except=self.minimize_except)
        self.params = {p:None for p in pars if not (p=="log_10(H0)" or p=='SNe Ia M_B')}
        if 'SNe Ia M_B' in pars:
          self.params["Mb"] = None

    def get_requirements(self):
        """
        Declare cosmological quantities needed from theory provider (CAMB/CLASS).
        """
        return {
            "H0": None,
            "angular_diameter_distance": {"z": self.zinterp}
        }

    def logp(self, **params_values):
        """
        Compute and return the log-likelihood value.
        """
        # Fetch requested quantities from theory provider
        H0 = self.provider.get_param("H0")

        # Personal grievance:
        # Cobaya is written in a way where the get_angular_diameter distance always expects to be executed over the same z array, which is a WILD assumption, TBH
        # It forces us to evaluate the function -- which is an interpolator already -- over a given fixed array, just to build an interpolator to pass on
        c = 299792.458
        cz_func = CubicSpline(self.zinterp,self.provider.get_angular_diameter_distance(self.zinterp)/self.zinterp*H0/c)
        kz_func = CubicSpline(self.zinterp,self.provider.get_angular_diameter_distance(self.zinterp)*(1+self.zinterp)**2/self.zinterp*H0/c)

        # Pass observables and sampling parameters to your custom evaluation function
        self.pip.adjust_equations({"cz_function":cz_func, "kz_function":kz_func})

        # This is only possible (and fast) because the covariance does NOT depend on cosmology. If this is not the case in your model, you HAVE to re-derive the icov (using e.g. self.pip.solve)
        pars, lkl, bf, cov = self.pip.get_likelihood(inv_covar = self.icov, minimize_except = self.minimize_except)

        param_vec = np.empty_like(pars, dtype=float)
        for ipar, par in enumerate(pars):
          if par == "log_10(H0)":
            param_vec[ipar] = np.log10(H0)
          elif par == 'SNe Ia M_B':
            param_vec[ipar] = params_values['Mb']
          else:
            param_vec[ipar] = params_values[par]

        loglike = lkl(param_vec)

        return loglike

