function [ITC_data, ERSP_data, freqs, nchan, frames, ntrial, erspboot_mat, ...
    itcboot_mat, p_ERSP, p_ITC, corr_mat, time_bins, bin_size, nbins, ...
    filenames] = hilbert_itc(data, srate, freq_range, freq_int, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: hilbert_itc
%       By: Kim Stachenfeld, Aug 9 2012                                         
%  Purpose: get phase and amplitude information in the time and frequency
%           domain of hilbert transformed data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Usage:
% >> [f_data, f_ITC, f_bootstrap] = plot_ITC(raw_data, channel, ...
%    freqs, times, ITC_data, boot_mat, p_mat, circ_corr_mat, time_bins);
% 
% >> [f_data, f_ITC, f_bootstrap] = plot_ITC(raw_data, channel, ...
%    freqs, times, ITC_data, boot_mat, p_mat, circ_corr_mat, time_bins ...
%    'key1', value1);
%
% Outputs:
%   ITC_data     = [nchan x nfreq x frames] complex values representing the
%                  intertrial coherence calculated from bandpass filtered, 
%                  Hilbert transformed data. ITC is calculated using the formula:
%
%                    ITC = 1/ntrial * sum(H_n(f,t)/|H_n(f,t)|) for n = 1:ntrial
%
%                  where H_n(f,t) represents the Hilbert transform of the data
%                  bandpassed for frequency f at time t over trial n of
%                  ntrial.
%   ERSP_data    = [nchan x nfreq x frames] 
%
%                   ERSP = 1/ntrial * sum(|H_n(f,t)|^2) for n = 1:ntrial
%
%   freqs        = [nfreq+1] frequencies at the boundary of each frequency bin
%                 (in Hz)
%   nchan        = number of channels
%   frames       = total number of time points
%   ntrial       = number of trials
%
% Optional Outputs (empty if optional inputs do not produce these outputs)
%   erspboot_mat = [nchan x nfreq x frames] std deviation of ERSP over 
%                  bootstrap distribution
%   itcboot_mat  = [nchan x nfreq x frames] std deviation of ERSP over 
%                  distribution made by resampling
%   p_ERSP       = [nchan x nfreq x frames]p-value of ERSP based on
%                  bootstrap distribution
%   p_ITC        = [nchan x nfreq x frames] p-value of ITC based on
%                  bootstrap distribution
%   corr_mat     = [nchan x nfreq x n_timebins] values from [0 - 1] 
%                  representing the average circular correlation, 
%
%                  Hilbert transformed data. ITC is calculated using the formula:
%
%                    corr = 1/ntrial^2 * sum(corr{ang[H_n(f,t)], ang[H_m(f,t)]}))
%                    for n = 1:ntrial, m = 1:ntrial, n~=m.
%                  where ang[H_n(f,t)] represents the Hilbert transform of the
%                  data for frequency f at time bin t over trial n of
%                  ntrial.
%
%   time_bins    = [n_timebins] fram
%   bin_size     = number of frames per time bin
%   nbins        = number of time bins
%   filenames    = filenames to which data is saved
%
% Required inputs:
% (Note - to omit any input and use default value, input empty matrix, [])
%   data           = [nchan x frames x ntrial] || [frames x ntrial] - Continuous
%                    data for ntrial trials over nchan channels
%   srate          = scalar||[] - sampling rate in Hz {default: 1000}
%   freq_range     = [minfreq maxfreq]
%   freq_int       = scalar < diff(freq_range) - size of intervals over which 
%                    data will be bandpass filter (for high frequency
%                    resolution, make freq_int small)
% Optional phase angle correlation inputs:
%   'angle_corr'   = [{1}||0] if 1, run phase angle correlation analysis
%   'n_timebins'   = integer - number of bins in the time domain for phase
%                    angle correlation analysis {100, unless 'timebin_size' 
%                    is specified}
%   'timebin_size' = size of bins in the time domain for phase angle 
%                    correlation analysis {floor(frames/timebin_size)} 
% Optional bootstrapping inputs:
%   'bootstrap'    = [{1}||0] if 1, run bootstrap analysis
%   'naccu'        = integer - number of accumulations {200}
%   'alpha'        = [0 0.5] - p-value rejection threshold (recommended: 
%                    alpha = 0.1 or alpha = 0.05) {0}
% Optional plotting parameters:
%   'plot_all'     = [{1}||0] if 1, plot all results
%   'plotrawdata'  = [{1}||0] if 1, plot raw data
%   'plotcalc'     = [{1}||0] if 1, plot ERSP, ITC, and phase angle 
%                    correlation calculation results
%   'plotboot'     = [{1}||0] if 1, plot bootstrap std deviation for ERSP
%                    and ITC
% Optional file save:
%   'save_figs'    = [1||{0}] - if 1, save figures
%   'save_data'    = [1||{0}] - if 1, save data
%   'save_all'     = [1||{0}] - if 1, save figures and data
%   'save_to_path' = string - filepath to which results will be saved {'.'}

g = finputcheck(varargin, ... 
  ...%PHASE ANGLE CORRELATION PARAMETERS
    {'angle_corr'    'real'     []                1;    ...
     'n_timebins'    'integer'  []                100;  ...
     'timebin_size'  'integer'  []                [];   ...
     %BOOTSTRAPPING PARAMETERS
     'bootstrap'     'real'     []                1;    ...
     'naccu'         'integer'  [0 Inf]           200;  ...
     'alpha'         'real'     [0 0.5]           0;    ...
     %PLOTTING PARAMETERS
     'plot_all'      'real'     []                1;    ...
     'plotrawdata'   'real'     []                1;    ...
     'plotcalc'      'real'     []                1;    ...
     'plotboot'       'real'     []               1;   ...
     'plotchan'      'integer'  [1:size(data,1)]  1;    ...
     'norm_ersp'     'real'     []                0;    ...
     %SAVE FILE PARAMETERS
     'save_all'      'real'     []                0;    ...
     'save_to_path'  'string'   []                './'; ...
     'save_figs'     'real'     []                0;    ...
     'save_data'     'real'     []                0     ...
     });

if g.save_all || g.save_figs || g.save_data
    if g.save_to_path(end) ~= '/'
        g.save_to_path = strcat(g.save_to_path, '/');
    end
end

filenames = [];
freqs = freq_range(1):freq_int:freq_range(2);
nfreq = length(freqs)-1;

dim = length(size(data));
   
if dim == 2
    data = permute(data, [3 1 2]);
else if dim > 3
        error('error - data matrix cannot have more than 3 dimensions.');
    end
end

nchan  = size(data, 1);
frames = size(data, 2);
ntrial = size(data, 3);
all_hf_data = zeros(nchan, nfreq, frames, ntrial);

for ch = 1:nchan
    for f = 1:nfreq
        %filter parameters
        locutoff = freqs(f+1);
        hicutoff = freqs(f);
        filtorder = min(floor(frames/3), 3*fix(srate/locutoff));
        if mod(filtorder,2) == 1
            filtorder = filtorder-1;
        end
        fprintf('Calculating intertrial coherence for Channel %0g, frequency %0g.\n', ch,f)
        for tr = 1:ntrial
            % bandpass filter data
            filt_data = eegfilt(squeeze(data(ch,:,tr)), srate, locutoff, 0, 0, filtorder);
            filt_data = eegfilt(filt_data, srate, 0, hicutoff, 0, filtorder);
            all_hf_data(ch, f, :, tr) = hilbert(filt_data);
        end
    end
end

%ERSP (power spectral density)
ERSP_data = mean(abs(all_hf_data).^2, 4);
for ch = 1:size(ERSP_data,1)
    for f = 1:size(ERSP_data,2)
        %ERSP_data(ch,f,:) = squeeze(ERSP_data(ch,f,:))/norm(squeeze(ERSP_data(ch,f,:)))^2;
    end
end



%Inter trial coherence
ITC_data = mean(all_hf_data./abs(all_hf_data), 4);


% Phase angle correlation
corr_mat = []; time_bins = 0; bin_size = 0; nbins = 0;
if g.angle_corr
    [corr_mat,time_bins,bin_size,nbins] = phase_angle_correlation(all_hf_data, ...
        nchan, frames, ntrial, nfreq, 'nbins', g.n_timebins, 'bin_size', g.timebin_size);
end

% Bootstrapping for ITC amplitude
erspboot_mat = []; itcboot_mat = []; p_ITC = []; p_ERSP = [];
if g.bootstrap
    [erspboot_mat, itcboot_mat, p_ITC, p_ERSP] = bootstrap_analysis...
        (all_hf_data, ERSP_data, ITC_data, g.naccu, nchan, nfreq, ntrial, frames);
end

% Plotting
if g.plot_all
    [f_raw, f_calc, f_boot] = plot_results(data, ERSP_data, ITC_data, ...
        erspboot_mat, itcboot_mat, p_ITC, p_ERSP, corr_mat, ...
        g.plotchan, freqs, -2000:2000, time_bins, 'alpha', g.alpha, ...
        'norm_ERSP', g.norm_ersp);
else
	 [f_raw, f_calc, f_boot] = plot_results(data, ERSP_data, ITC_data, ...
         erspboot_mat, itcboot_mat, p_ITC, p_ERSP, corr_mat, g.plotchan,...
         freqs, t, time_bins,'plot_raw', g.plotrawdata, 'plot_calc', g.plotcalc, ...
        'plot_boot', plotboot, 'norm_ERSP', g.norm_ersp, 'alpha', g.alpha);
end

% if data = [frames x ntrial], inverse permute results to 2D
if dim == 2
    ERSP_data = ipermute(ERSP_data, [3 1 2]);
    ITC_data = ipermute(ITC_data, [3 1 2]);
    corr_mat = ipermute(corr_mat, [3 1 2]);
    erspboot_mat = ipermute(erspboot_mat, [3 1 2]);
    itcboot_mat = ipermute(itcboot_mat, [3 1 2]);
    p_ITC = ipermute(p_ITC, [3 1 2]);
    p_ERSP = ipermute(p_ERSP, [3 1 2]);
end

% Save
if g.save_all
    g.save_figs = 1;
    g.save_data = 1;
end

if g.save_data
    filenames = {'ITC_data';'ERSP_data';'freqs';'nchan';'frames';'ntrial';
        'erspboot_mat';'itcboot_mat';'p_ERSP';'p_ITC';'corr_mat';'time_bins';
        'time_bins';'bin_size';'nbins'};
    for i = 1:length(filenames)
        save(strcat(g.save_to_path, filenames{i}), evalc(filenames{i}));
    end
end

if g.save_figs
    filenames = {filenames; 'rawdata_fig'; 'itc_ersp_pac_fig';'bootstrap_fig'};
    print(f_raw, strcat(g.save_to_path, 'rawdata_fig'), '-dpng');
    print(f_calc, strcat(g.save_to_path, 'itc_ersp_pac_fig'), '-dpng');
    print(f_boot, strcat(g.save_to_path, 'bootstrap_fig'), '-dpng');
end



fprintf('\nCONGRATULATIONS! You now know the ERSP, ITC, \nand maybe even the phase-angle correlation of your data.\n')
