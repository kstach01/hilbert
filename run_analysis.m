clear;
clc;
close all;

load ex_KI2;

%Set filtering parameters
freq_range = [1 51];
freq_int = 0.5;
alpha = 0.1;
naccu = 200;
channel = 1; %perform analysis on first channel only for demonstration purposes

%to call hilbert_itc, 
[ITC_data, ERSP_data, freqs, nchan, frames, ntrial, erspboot_mat, ...
    itcboot_mat, p_ERSP, p_ITC, corr_mat, time_bins, bin_size, nbins] ...
    = hilbert_itc(data(channel,:,:), srate, freq_range, freq_int, 'alpha', ...
    alpha, 'angle_corr',0,'bootstrap',1);

% for information on functions and optional inputs, enter into the command line:
% >> help hilbert_itc
