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

related works: 

“Maximum Entropy Method without False Peaks with Exact Numerical Equation”, J. Phys.: Conf. Ser., vol. 1438, 012031 (6pp), 2020
https://iopscience.iop.org/article/10.1088/1742-6596/1438/1/012031

"Refreshing idea on Fourier analysis," Proc. 21st IEEE CSPA, pp.1-4, 2025
https://arxiv.org/abs/2501.03514
https://doi.org/10.1109/CSPA64953.2025.10933386

QR for this page
![qrcode_github com](https://github.com/user-attachments/assets/c9e5970a-99ba-45f9-8667-b3b2f394adc8)
