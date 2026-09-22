"""
Simple standalone test / example usage script for the H0DN_like likelihood.
Run from anywhere after `pip install -e .`:
    python H0DN_like/example.py
"""
if __name__ == "__main__":

  from cobaya.run import run

  info = {
      "likelihood": {
          "H0DN_like": {
              #"modify_config": {"sn1a_hf_file" : "", "sn1a_hf_cov_file":"", "alpha_sn1a":0.714, "alpha_sn1a_error":1000.0},
          },
          #"sn.pantheonplus" : {'use_abs_mag':True,
          #                     #'zmax':0.15
          #                     },
      },
      "theory": {
          "classy": {
              "ignore_obsolete":True
          }
      },
      "params": {
          # Sampled cosmological parameters
          "H0": {"prior": {"min": 50, "max": 90}, "proposal": 1.0,'ref':73.5, "latex": r"H_0"},
          "omega_b": 0.02233,
          #"Omega_m": {"prior": {"min": 0.2, "max": 0.4}, 'ref':0.35, "proposal": 0.005, "latex": r"\Omega_m"},
          "Omega_m":0.3,
          "Mb" : {"prior": {"min":-20, "max":-19}, 'ref':-19.253, "proposal": 0.01, "latex": r"M_B"},
          "As": 2.1e-9,
          "ns": 0.965,
          "tau": 0.054,
      },
      "sampler": {
          "mcmc": {
              #"max_samples": 1000,
              "max_samples": 1_000_000,
              "Rminus1_stop": 0.005,
          }
      },
      "output": "chains/test_h0_run3",
      "resume": True,
  }

  updated_info, sampler = run(info)
  
if __name__ == "__main__":
  main()
