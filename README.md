Hilbert Transform for Brain Wave Analysis
=======

__Author:__ Kimberly Stachenfeld  
__Email:__ kim.stachenfeld@gmail.com

##Purpose
Analyze EEG or LFP data via Hilbert Transform

<!---
##Data Preprocessing
What form should the data be in?
-->

##Contents

`run_analysis`: Example on how to call `hilbert_stats`
`ex_KI2.mat`: Example dataset
`hilbert_stats`: Returns ITC, ERSP, and phase angle correlation and relevant bootstrap statistics and plots from hilbert transform data bandpassed at frequency intervals.
`bootstrap_analysis`: Calculates bootstrap and p-value by resampling for ITC and ERSP
`phase_angle_correlation`: Calculates average correlation among phase angles for different trials
`plot_results`: plot raw data, ERSP, ITC, phase angle correlation, and bootstrap results

##Dependencies
###*circstat* toolbox  
Philipp Behrens. [A Matlab Toolbox for Circular Statistics](http://www.jstatsoft.org/v31/i10). _Journal of Statistical Software_, 2009, __31__(10).
####Required from circStat package:  
- `circ_r`: Resultant vector length  
- `circ_mean`: Mean direction of a sample of circular data  
- `circ_confmean`: Confidence intervals for mean direction  
- `circ_corrcc`: Circular-circular correlation coefficient  

###*EEGLAB* package
Arnaud Delorme, Scott Makeig. [EEGLAB: an open source toolbox for analysis of single-trial EEG dynamics](http://sccn.ucsd.edu/~scott/pdf/EEGLAB04.pdf). *Journal of Neuroscience Methods*, 2004, **134**:9-21. 
####Required from eeglab package:
- `eegfilt`: Low/High pass filter data
<!--- note to kim: sub out with function to bandpass filter without changing phase information -->


