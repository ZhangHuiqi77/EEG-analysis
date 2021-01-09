function []= spectrogramFunc(x,dt)
% huiqi
% 02/01/2021
% this function generate a power spectrogram with colorbar(jet mode)
% Input:
%   x:data, should be a row vector
%   dt: sampling interval(s)
% Output:
%   none values
    Fs = 1/dt;					%Define the sampling frequency.
    interval = round(Fs);		%Specify the interval size.
    overlap = round(Fs*0.95);	%Specify the overlap of intervals.
    nfft = round(Fs);			%Specify the FFT length.

    %Compute the spectrogram,
    [S,F,T,P] = spectrogram(x-mean(x),interval,overlap,nfft,Fs);
    figure
    % imagesc(T,F,10*log10(P))	%... and plot it,
    imagesc(T,F,P)
    colorbar					%... with a colorbar,
    axis xy						%... and origin in lower left, 
end
