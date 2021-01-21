function [lags,crosscorr,max_crosscorr_lag] = amp_crosscorr_trials(x,y,lags_N,Fs,Freq1,Freq2)
% huiqi
% 20/01/2021
% amp_crosscorr_trials is for compute the amplitude crosscorrelation for a 
% series of epoches, which may indicate the direction of information flow. 
% called function: amp_crosscorr
% USAGE: [lags,crosscorr,max_crosscorr_lag,ind] = amp_crosscorr_trials(x,y,lags_N)
% Input:
%   x: signal, should be a vector.
%   y: signal, should be the same size of x.
%   lags_N: the limitation of lags[ms], usually use 100ms.
%   Fs: sampling frequency
%   Freq1: frequency limitation-low
%   Freq2: frequency limitation-high, must bigger than Freq1
%   (Freq REF: delta(1-4),
%              theta(7-12),
%              beta(13-30),
%              low gamma(30-50),
%              high gamma(50-100)).
% Output:
%   lags: lag vector for each trial, actually the same for all trials.
%   crosscorr: transfer the result calculated by amp_crosscorr function.
%   max_crosscorr_lag
%   ind: the index of sorted crosscorr
    ntrials = size(x,1);
    lags = zeros(ntrials,2*lags_N+1);
    crosscorr = zeros(ntrials,2*lags_N+1);
    max_crosscorr_lag = zeros(ntrials,1);
    for k = 1:ntrials
        [lags(k,:),crosscorr(k,:),max_crosscorr_lag(k)] = amp_crosscorr(x(k,:),y(k,:),Fs,Freq1,Freq2);
    end

    figure;
    [max_crosscorr_lag,ind] = sort(max_crosscorr_lag);
    crosscorr = crosscorr(ind,:);
end