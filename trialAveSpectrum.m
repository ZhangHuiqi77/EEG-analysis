function [faxis, Sxxm] = trialAveSpectrum(E1,dt)
% huiqi
% 02/01/2021
% This function calculate the average spectrum of a series trials. 
% Input
%   E1: data structure should be every row represent a trial
%   dt represent the sampling interval(s)
% Output
%   faxis: frequency axis
%   Sxxm: mean spectrum of all trials, every trial is calculated based on
%   fft and autocorrelation.

    K = size(E1,1);				%Define the number of trials.
    N = size(E1,2);				%Define the number of time indices.
    t = (dt:dt:N*dt);           %Generate time viarable
    T  = t(end);				%Define the duration of data.

    Sxx = zeros(K,N);		%Create variable to store each spectrum.
    for k=1:K					%For each trial,
        x = E1(k,:);			%... get the data,
        xf  = fft(x-mean(x));	%... compute Fourier transform,
        Sxx(k,:) = 2*dt^2/T *(xf.*conj(xf));%... compute spectrum.
    end
    Sxx = Sxx(:,1:N/2+1);		%Ignore negative frequencies.
    Sxxm = mean(Sxx,1);			%Average spectra over trials.

    df = 1/max(T);				%Define frequency resolution,
    fNQ = 1/dt/2;				%... and Nyquist frequency.
    faxis = (0:df:fNQ);			%... to construct frequency axis.

end