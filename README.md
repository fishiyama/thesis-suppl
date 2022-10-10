# thesis-suppl

### supplemental code [tutorial.jl](https://github.com/fishiyama/thesis-suppl/blob/main/tutorial.jl) for my thesis https://tsukuba.repo.nii.ac.jp/records/2000757.

we demonstrate our nonlinear method of time-frequency analysis with a Julia code.

our method is a kind of mode decomposition with general complex functions.

our Julia code works as follows:

time series for analysis. (12 samples)

![Fig1 timeseries for anaysis](https://user-images.githubusercontent.com/111185366/192537801-fb49c8c7-c94a-47d9-b027-faaf73d3d53d.png)

obtained results.

![Fig2 obtained results](https://user-images.githubusercontent.com/111185366/192537829-617fd562-6fa0-4766-80ac-1178ae0e1245.png)

no.1 backgound noise, no.2&3 frequency of sinusoidal term frq=0.1, no.4 constant term amp=0.01
 
you will find following values with runnning our code.
![image](https://user-images.githubusercontent.com/111185366/193454046-6604d2aa-f803-40c6-88f4-d7c13eb85da3.png)
please remind the "theoretical limit of time-frequency resolution", and compare it with this code.
which do you think is the reality?
