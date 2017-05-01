# GMM_IC_DP

"Initializing EM algorithm for univariate Gaussian, multi-component, heteroscedastic mixture models by dynamic programming partitions"

by Andrzej Polanski, Michal Marczyk, Monika Pietrowska, Piotr Widlak, Joanna Polanska

SUPPLEMENTARY MATERIALS:
Matlab scripts and function for performing comparisons of partitioning
algorithms E-Q, H-clu-c, H-clu-a, DP-Q4 for the data described as
Group 4 in section 6.2. Computations are started by launching Matlab
script 
                partitions_em_demo

One stochastic simulation experiment is done (including three steps
1-3 listed in section 6.2). Results of computations are shown by 
plots of partitions (figure 1) and data histograms versus estimated 
probability density functions (figure 2). Values of errors 
and likelihoods are also reported. By modifications of the Matlab code 
other computational scenarios (for simulated data) can be also realized.


LIST OF FILES:
partitions_em_demo: script file for demonstratig comparison of algorithms 
                    for comparisons of initialization methods for EM algorithm
g_mix_gen: function, generates mixture sample
h_clu_a: function, sample partition by average linkage hierarchical clustering
h_clu_c: function, sample partition by complete linkage hierarchical clustering
dyn_pr_split: function, computes sample partition by dynamic programming  
g_mix_est_fast_lik: function, performs EM iterations 
draw_part: function, draws partition
draw_hist_pdf: function, draws data histogram versus pdf of the estimated mixture
comp_errors: function for computing error (scaled difference between true and 
             estimated parameters).


USAGE:
Copy all files to one folder, launch: partitions_em_demo

[Andrzej Polanski; andrzej.polanski@polsl.pl][Michał Marczyk; michal.marczyk@polsl.pl]# GMM_IC_DP
