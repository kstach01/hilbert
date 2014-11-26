function [corr, time, bin_size, nbins] = phase_angle_correlation ...
    (hf_data, nchan, frames, ntrial, nfreq, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function: phase_angle_correlation(hf_data,nchan,frames,ntrial,nfreq,...)
%      By: Kim Stachenfeld, Aug 6 2012                                         
% Purpose: Divide complex data into time bins and find the average 
%          circular Pearson Correlation Coefficient among trials for each 
%          time bin.
%Requires: CircStat package
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Usage:
% >> [circ_corr,p_value,time,bin_size,nbins] = phase_angle_correlation ...
%            (hf_data, nchan, frames, ntrial, nfreq);
%
% >> [circ_corr p_value time bin_size nbins] = phase_angle_correlation ...
%            (hf_data, nchan, frames, ntrial, nfreq, ..., 'key1', value1);
%
% Outputs:
%   circ_corr = [nchan x nfreq x nbins] Average circular correlation 
%                 coefficient for each channel, frequency bin, and time bin
%   time      = vector of frames at which each bin starts.
%   bin_size  = size of each bin
%   nbins     = number of bins
%
% Required inputs:
%   hf_data   = [nchan x nfreq x frames x ntrial] matrix containing
%               complex data for filtered and hilbert transformed data 
%               at all channels, frequency bins, times, and trials
%   nchan     = number of channels {default: size(hf_data, 1)}
%   frames    = number of frames/time points {default: size(hf_data, 3)}
%   ntrial    = number of trials {default: size(hf_data, 4)}
%   nfreq     = number of frequency bins {default: size(hf_data, 2)}
%
% Optional inputs:  Value                                 {default}
%   nbins     = number of bins into which to divide frames. A large number
%   will provide high time resolution but poor correlation resolution.
%   bin_size  = size of bins
%
% See also: hilbert_itc, circ_corrcc

g = finputcheck(varargin, ...
    {'bin_size'  'integer'   [0 frames]   []; ...
     'nbins'     'integer'   [0 frames]   []
     });

% Check for empty variables
if isempty(nchan)
    nchan = size(hf_data, 1);
end
if isempty(nfreq)
    nfreq = size(hf_data, 2);
end
if isempty(frames)
    frames = size(hf_data, 3);
end
if isempty(ntrial)
    ntrial = size(hf_data, 4);
end 
if isempty(g.nbins)
    if ~isempty(g.bin_size)
        g.nbins = floor(frames/g.bin_size);
    else
        g.nbins = floor(frames/10);
    end
end
if isempty(g.bin_size)
    g.bin_size = floor(frames/g.nbins);
end
% if ~matlabpool('size')
%     matlabpool open;
% end
%Check that bin size and bin number makes sense
if g.bin_size*g.nbins > frames
    error('bin_size * number_of_bins is greater than the number of frames');
end

den = sum(1:ntrial-1);
hf_data = angle(hf_data);
% Allocate memory for results
corr = zeros(nchan, nfreq, g.nbins);
% p_val = zeros(nchan, nfreq, g.nbins);
bin_data = zeros(g.bin_size, ntrial);
bins = reshape([g.bin_size*(0:g.nbins-1)+1; g.bin_size*(1:g.nbins)], g.nbins*2, 1);

fprintf('Calculating correlation of phase angle among trials.\n');

for ch = 1:nchan
    for f = 1:nfreq
        fprintf('   Analyzing channel %0g, frequency %0g.\n', ch, f);
        fprintf('\t Performing correlation for bin (of %0g): [0', g.nbins);
        N = 1;
        for b = 1:2:2*g.nbins
            a = b+1;
            if ~mod(N, floor(g.nbins/5))
                fprintf(' %0g', N); %print progress updates
            end
            bin_data(:,:) = squeeze(hf_data(ch, f, bins(b):bins(a), :));
            %average circular correlation and p_value among trials
            for tr1 = 1:ntrial-1
                for tr2 = tr1+1:ntrial
                	c = circ_corrcc(bin_data(:,tr1), bin_data(:,tr2));
                    corr(ch, f, a/2) = corr(ch, f, a/2) + c;
                    %p_val(ch, f,a/2) = p_val(ch, f, a/2) + p;
                end
            end
            N = N+1;
        end
        fprintf(']\n');
    end
end

%Outputs
corr = corr/den;
time = bins(1:2:end);
bin_size = g.bin_size;
nbins = g.nbins;
end
