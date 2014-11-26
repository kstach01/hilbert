function [erspboot_mat, itcboot_mat, p_ITC, p_ERSP] = bootstrap_analysis ...
    (hf_data, ERSP_data, ITC_data, naccu, nchan, nfreq, ntrial, frames)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BOOTSTRAP_ANALYSIS
% BY: Kim Stachenfeld
% UPDATED: August 7, 2012
% PURPOSE: perform resampling to estimate ITC bootstrap and p-value.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage:
% >> % To use all default parameters:
% >> [erspboot_mat, itcboot_mat, p_ersp, p_ITC] = bootstrap_analysis(hf_data);
%
% >> % To manually set some parameters and use other defaults:
% >> [erspboot_mat, itcboot_mat, p_ersp, p_ITC] = bootstrap_analysis ...
% (hf_data, ITC_data, [], nchan);
%
% >> % To manually set all parameters:
% >> [erspboot_mat, itcboot_mat, p_ersp, p_ITC] = bootstrap_analysis ...
%            (hf_data, ITC_data, naccu, nchan, nfreq, ntrial, frames);
%
% Outputs:
%      erspboot_mat  = [nchan x nfreq x frames] matrix containing bootstrap
%                      for each channel, frequency bin, and frame.
%      itcboot_mat   = [nchan x nfreq x frames] matrix containing bootstrap
%                      for each channel, frequency bin, and frame.
%      p_ERSP        = [nchan x nfreq x frames] matrix containing p-value for
%                      each channel, frequency bin, and frame. The p-value is 
%                      the probability of finding a value as extreme as that 
%                      found.
%      p_ITC         = [nchan x nfreq x frames] matrix containing p-value for
%                      the ITC at each channel, frequency bin, and frame.
%
% Required Inputs:
% Note: to use default value, enter [] or omit entry
%      HF_DATA   = [nchan x nfreq x frames x ntrial] matrix containing 
%                  hilbert transformed data that was bandpass filtered 
%                  into nfreq frequency bins{default: error}       
%      ERSP_data = 
%      ITC_data  = [nchan x nfreq x frames x ntrial]{default: calculate ITC}
%      naccu     = number of accumulations in resampling. Note that 
%                  large values will produce less variability, but will 
%                  significantly increase run time {def: 200}
%      nchan     = number of channels (scalar) {def: size(hf_data,1)
%      nfreq     = number of frequency bins (scalar)
%                  {default: size(hf_data,2)}
%      ntrial    = number of trials (scalar) {def: size(hf_data,4)}
%      frames    = number of frames in epoch (scalar) 
%                  {default: size(hf_data, 3)}
%
% See also: hilbert_itc

%Set unfilled inputs to default values
if isempty(hf_data)
    error('raw data matrix cannot be empty');
end
if isempty(ITC_data)
    ITC_data = abs(mean(hf_data./abs(hf_data), 4));
else if ~isempty(find(imag(ITC_data) ~= 0, 1))
        ITC_data = abs(ITC_data);
    end
end
if isempty(nchan)
    nchan = size(hf_data,1);
end
if isempty(nfreq)
    nfreq = size(hf_data,2);
end
if isempty(frames)
    frames = size(hf_data,3);
end
if isempty(ntrial)
    ntrial = size(hf_data,4);
end
if isempty(naccu)
    naccu = 200;
end

%bootstrapping data structures
erspboot_mat = zeros(nchan, nfreq, frames);
itcboot_mat = zeros(nchan, nfreq, frames);

p_ITC = zeros(nchan, nfreq, frames);
p_ERSP = zeros(nchan, nfreq, frames);

temp_ITC = zeros(frames, naccu);
tempERSP = zeros(frames, naccu);

t = randi(ntrial, ntrial, 1);

% Bootstrap time!
fprintf('Performing bootstrap resampling with %0g accumulations.\n', naccu);
for ch = 1:nchan
    for f = 1:nfreq
        fprintf('   Analyzing channel %0g, frequency %0g.\n', ch, f);
        fprintf('\t Performing accumulation (of %0g): [0', naccu);
        n = 1;
        %Perform accumulations
        while n <= naccu
            if ~mod(n, floor(naccu/5))
                fprintf(' %0g', n); %print progress updates
            end
            temp_ITC(:,n) = squeeze(mean( hf_data(ch,f,:,t)./abs(hf_data(ch,f,:,t)), 4 ));
            tempERSP(:,n) = squeeze(mean( mean(abs(hf_data(ch,f,:,t)).^2, 4), 4 ) );
            t = randi(ntrial, ntrial, 1);
            n = n+1;
        end
        fprintf(']\n');
        % P-value: probability of finding a value as extreme as that found.
        % It will be the min( fraction of values that are larger, fraction
        % of values that are smaller )
        p_ITC(ch, f, :)  = min( mean( bsxfun( @ge, abs(temp_ITC), squeeze(ITC_data(ch,f,:)) ), 2), ...
                                mean( bsxfun( @le, abs(temp_ITC), squeeze(ITC_data(ch,f,:)) ), 2) );
        p_ERSP(ch, f, :) = min( mean( bsxfun( @ge, abs(tempERSP), squeeze(ERSP_data(ch,f,:)) ), 2), ...
                                mean( bsxfun( @le, abs(tempERSP), squeeze(ERSP_data(ch,f,:)) ), 2) );
        % Bootstrap: standard deviation of distribution produced by
        % resampling
        itcboot_mat(ch, f, :) = std(temp_ITC, [], 2);
        erspboot_mat(ch, f, :) = std(tempERSP, [], 2);
    end
end



