# CF-LPZF

# Local Partial Zero-Forcing Combining for Cell-Free Massive MIMO Systems

This is a code package is related to the following scientific article:

Jiayi Zhang, Jing Zhang, Emil Björnson, ,and Bo Ai, "[Local Partial Zero-Forcing Combining for Cell-Free Massive MIMO Systems](https://ieeexplore.ieee.org/document/9529197)," *IEEE Transactions on Communications*, vol. 69, no. 12, pp. 8459-8473, December 2021.

The package contains a simulation environment, based on Matlab, that reproduces some of the numerical results and figures in the article. *We encourage you to also perform reproducible research!*


## Abstract of Article

Cell-free massive multiple-input multiple-output (MIMO) provides more uniform spectral efficiency (SE) for users (UEs) than cellular technology. The main challenge to achieve the benefits of cell-free massive MIMO is to realize signal processing in a scalable way. In this paper, we consider scalable full-pilot zero-forcing (FZF), partial FZF (PFZF), protective weak PFZF (PWPFZF), and local regularized ZF (LRZF) combining by exploiting channel statistics. We derive closed-form expressions of the uplink SE for FZF, PFZF, and PWPFZF combining with large-scale fading decoding over independent Rayleigh fading channels, taking channel estimation errors and pilot contamination into account. Moreover, we investigate the impact of the number of pilot sequences, antennas per AP, and APs on the performance. Numerical results show that LRZF provides the highest SE. However, PWPFZF is preferable when the number of pilot sequences is large and the number of antennas per AP is small. The reason is t at PWPFZF has lower computational complexity and the SE expression can be computed in closedform. Furthermore, we investigate the performance of PWPFZF combining with fractional power control and the numerical results show that it improves the performance of weak UEs and realizes uniformly good service for all UEs in a scalable fashion.

## Content of Code Package

The package generates the simulation SE results which are used in Figure 1, Figure 2, Figure 3, Figure 4, Figure 5, and Figure 6. To be specific:

- `mRZF`: Compute simulation SE results while using mRZF;
- `simu_anal_MR_FZF_LPZF`: Compute the SE by Monte Carlo simulation and analytical expression using MR, FZF, and LPZF;
   - `functionComputeSE_AP_uplink_anal_ADC_FZF`
   - `functionComputeSE_AP_uplink_anal_ADC_LPZF_final`
   - `functionComputeSE_AP_uplink_anal_ADC_LPPZF`
   - `functionComputeSE_AP_uplink_simu_ADC_FZF`
   - `functionComputeSE_AP_uplink_simu_ADC_LPZF`
   - `functionComputeSE_AP_uplink_simu_ADC_LPPZF`
- `generateSetup_new`: Generate spatial correlation matrix and large-scale fading coefficients;
  - `functionRlocalscattering`: Generate the spatial correlation matrix for the local scattering model;
- `MR_MMSE_FZF_LPZF_LPPWZF_LPPSZF_simu_p_equal_K10': Compute simulation SE results while using MMSE, FZF, LPZF, LPPWZF, and LPPSZF;
- `functionUEgrouping`: Perform UE grouping;
- `fractionalPowerControl`: Perform uplink fractional power control; 
- `functionChannelEstimates`: Perform MMSE channel estimation;
- `functionComputeSE_AP_uplink_lowADC_FZF_RZF_final`: Compute simulation SE results while using FZF and RZF;
- `functionComputeSE_AP_uplink_lowADC_LPZF_LPPZF_final`: Compute simulation SE results while using LPZF and LPPZF;

See each file for further documentation.


## License and Referencing

This code package is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our original article listed above.
