function [faxis,cohr] = cohrFunc(E1,E2,Fs)
% huiqi
% 03/01/2021
% This function is for computing multi-trial average coherence between two electrodes signals.
% Also works for single trial coherence.

% Input:
%   E1: The first electrode signal. Each row is a trial.
%   E2: The second electrode signal. Shoud be the same size with E1.
%   Fs: sampling frequency.

% Output:
%   faxis: frequency axis.
%   cohr: computed cohr values.

K = size(E1,1);			%Define the number of trials.
N = size(E1,2);			%Define the number of indices per trial.
dt = 1/Fs;              %Define the sampling interval.
t = (dt:dt:N*dt);       %Define the timeline according to epoch.
T  = t(end); 			%Define the duration of data.

Sxx = zeros(K,N);		%Create variables to save the spectra,
Syy = zeros(K,N);
Sxy = zeros(K,N);
for k=1:K				%... and compute spectra for each trial.
    x=E1(k,:)-mean(E1(k,:));    % offset bias of E1.
    y=E2(k,:)-mean(E2(k,:));    % offset bias of E2.
    Sxx(k,:) = 2*dt^2/T * (fft(x) .* conj(fft(x)));
    Syy(k,:) = 2*dt^2/T * (fft(y) .* conj(fft(y)));
    Sxy(k,:) = 2*dt^2/T * (fft(x) .* conj(fft(y)));
end

Sxx = Sxx(:,1:N/2+1);	%Ignore negative frequencies.
Syy = Syy(:,1:N/2+1);
Sxy = Sxy(:,1:N/2+1);

Sxx = mean(Sxx,1);		%Average the spectra across trials.
Syy = mean(Syy,1);
Sxy = mean(Sxy,1);		%... and compute the coherence.
cohr = abs(Sxy) ./ (sqrt(Sxx) .* sqrt(Syy));

df = 1/max(T);			%Determine the frequency resolution.
fNQ = 1/dt/2;			%Determine the Nyquist frequency,
faxis = (0:df:fNQ);		%... and construct frequency axis.
