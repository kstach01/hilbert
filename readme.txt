hilbertEEG
=======================

Purpose: Analyze EEG or LGP data via Hilbert Transform

 Author: Kimberly Stachenfeld
  Email: kim.stachenfeld@gmail.com

Contents:

run_analysis                Example on how to call hilbert_stats
ex_KI2.mat                  Example dataset
hilbert_stats               Returns ITC, ERSP, and phase angle correlation and relevant bootstrap statistics and plots
			    from hilbert transform data bandpassed at frequency intervals.
bootstrap_analysis          Calculates bootstrap and p-value by resampling for ITC and ERSP
phase_angle_correlation     Calculates average correlation among phase angles for different trials
plot_results                plot raw data, ERSP, ITC, phase angle correlation, and bootstrap results

Required from circStat package:
circ_r 		            Resultant vector length

circ_mean 		    Mean direction of a sample of circular data

circ_confmean 		    Confidence intervals for mean direction


circ_corrcc		    Circular-circular correlation coefficient

 
A Matlab Toolbox for Circular Statistics, Journal of Statistical Software, Volume 31, Issue 10, 2009.
	http://www.jstatsoft.org/v31/i10



Required from eeglab package:
eegfilt                     Low/High pass filter data

Delorme, A., Makeig, S. (in press) EEGLAB: An open source toolbox for analysis of single-trial 
EEG dynamics including independent component analysis. Journal of Neuroscience Methods.
	http://sccn.ucsd.edu/eeglab/download/eeglab_jnm03.pdf