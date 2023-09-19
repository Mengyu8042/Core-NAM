# Nonparametric Additive Models for Billion Observations
This repository includes the implementation of our work **"Nonparametric Additive Models for Billion Observations"**.


## Introduction
A brief introduction about the folders and files:
* `data/`: path for storing downloaded data;
* `source/`: the source file for sparse matrix operations;
* `utility.R`: containing useful functions used in the following scripts;
* `simu_vary_r.R`: simulation script; comparison of FULL, UNIF, LowCon, and CORE w.r.t. MSE and PMSE versus increasing subsample size $r$;
* `simu_vary_n.R`: simulation script; comparison of FULL, UNIF, LowCon, and CORE w.r.t. MSE and PMSE versus increasing full sample size $n$;
* `simu_sensitivity.R`: simulation script; sensitivity analysis by varying the signal-to-noise ratio and basis dimension;
* `simu_ineff.R`: simulation script; empirically verifying the asymptomatic optimality of CORE-GCV;
* `real_data_example.R`: real data analysis script; comparison of FULL, UNIF, LowCon, and CORE w.r.t. PMSE and PMAE versus increasing subsample size $r$.


## Reproducibility
For simulation studies in Section 5,
* you can run `simu_vary_r.R` to reproduce the results in Figure 3;
* you can run `simu_vary_n.R` to reproduce the results in Figure 4;
* you can run `simu_sensitivity.R` to reproduce the results in Figures 6 and 7;
* you can run `simu_ineff.R` to reproduce the results in Figure 8;

For the real data example in Section 6,
* first, download daily TCO data (.nc files as follows) during 1978--2019 from the link https://zenodo.org/record/7447660#.ZCJWMspBzEY and store them in the `data/` path;

        "NIWA-BS_CombinedTCO_V3.4.1_1978_Daily_Unpatched.nc"
        "NIWA-BS_CombinedTCO_V3.4.1_1979_Daily_Unpatched.nc"
                                ...
        "NIWA-BS_CombinedTCO_V3.4.1_2019_Daily_Unpatched.nc"

* then, run `real_data_example.R` to reproduce the results in Figure 9.