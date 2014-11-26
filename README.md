Hilbert Transform for Brain Wave Analysis
=======

Purpose: Analyze EEG or LFP data via Hilbert Transform

 Author: Kimberly Stachenfeld
  Email: kim.stachenfeld@gmail.com

Contents:

run\_analysis                Example on how to call hilbert\_stats
ex\_KI2.mat                  Example dataset
hilbert\_stats               Returns ITC, ERSP, and phase angle correlation and relevant bootstrap statistics and plots
			    from hilbert transform data bandpassed at frequency intervals.
bootstrap\_analysis          Calculates bootstrap and p-value by resampling for ITC and ERSP
phase\_angle\_correlation     Calculates average correlation among phase angles for different trials
plot\_results                plot raw data, ERSP, ITC, phase angle correlation, and bootstrap results

Required from circStat package:
circ\_r 		            Resultant vector length

circ\_mean 		    Mean direction of a sample of circular data

circ\_confmean 		    Confidence intervals for mean direction


circ\_corrcc		    Circular-circular correlation coefficient

 
A Matlab Toolbox for Circular Statistics, Journal of Statistical Software, Volume 31, Issue 10, 2009.
	http://www.jstatsoft.org/v31/i10



Required from eeglab package:
eegfilt                     Low/High pass filter data
\# note to kim: sub out with function to bandpass filter without changing phase information

Delorme, A., Makeig, S. (in press) EEGLAB: An open source toolbox for analysis of single-trial 
EEG dynamics including independent component analysis. Journal of Neuroscience Methods.
	http://sccn.ucsd.edu/eeglab/download/eeglab_jnm03.pdf
