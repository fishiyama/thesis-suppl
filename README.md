# thesis-suppl

### supplemental code [tutorial.jl](https://github.com/fishiyama/thesis-suppl/blob/main/tutorial.jl) for my thesis https://tsukuba.repo.nii.ac.jp/records/2000757.

we demonstrate our nonlinear method of time-frequency analysis with a [Julia](https://julialang.org/) code.

our method is a kind of mode decomposition with general complex functions.

our Julia code works as follows:

time series for analysis. (12 samples)

![Fig1 timeseries for anaysis](https://user-images.githubusercontent.com/111185366/192537801-fb49c8c7-c94a-47d9-b027-faaf73d3d53d.png)

obtained results.

![Fig2 obtained results](https://user-images.githubusercontent.com/111185366/192537829-617fd562-6fa0-4766-80ac-1178ae0e1245.png)

no.1 backgound noise, no.2&3 frequency of sinusoidal term frq=0.1, no.4 constant term amp=0.01
 
you will find following values with runnning our code.
![image](https://user-images.githubusercontent.com/111185366/193454046-6604d2aa-f803-40c6-88f4-d7c13eb85da3.png)
please compare our results with the "theoretical limit of time-frequency resolution" of conventional methods.

### related works: 

#### "Refreshing idea on Fourier analysis," Proc. 21st IEEE CSPA, pp.1-4, 2025
#### https://arxiv.org/abs/2501.03514 / https://doi.org/10.1109/CSPA64953.2025.10933386
``The "theoretical limit of time-frequency resolution in Fourier analysis" is thought to originate in certain mathematical and/or physical limitations. This, however, is not true. The actual origin arises from the numerical (technical) method deployed to reduce computation time. ...''

#### “Maximum Entropy Method without False Peaks with Exact Numerical Equation”, J. Phys.: Conf. Ser., vol. 1438, 012031 (6pp), 2020
#### https://iopscience.iop.org/article/10.1088/1742-6596/1438/1/012031
``The standard numerical maximum entropy method (MEM) still uses the Yule-Walker equation which contains rough approximation by Walker. The commonly used numerical equation contains additional modifications to reduce calculation cost. Nowadays, we have powerful computers and there is no reason to use the modified equation. ...''

### applications:

#### "Instantaneous Spectra Analysis of Pulse Series - Application to Lung Sounds with Abnormalities," to appear Proc. 22nd IEEE CSPA.
#### https://arxiv.org/abs/2602.03680
``The origin of the "theoretical limit of time-frequency resolution of Fourier analysis" is from its numerical implementation, especially from an assumption of "Periodic Boundary Condition (PBC)," which was introduced a century ago. We previously proposed to replace this condition with "Linear eXtrapolation Condition (LXC)," which does not require periodicity. This feature makes instantaneous spectra analysis of pulse series available, which replaces the short time Fourier transform (STFT). We applied the instantaneous spectra analysis to two lung sounds with abnormalities (crackles and wheezing) and to a normal lung sound, as a demonstration. Among them, crackles contains a random pulse series. ...''


QR for this page
![qrcode_github com](https://github.com/user-attachments/assets/c9e5970a-99ba-45f9-8667-b3b2f394adc8)
