function [seg]= segdata(x,t0,Fs,dur)
% huiqi
% 03/01/2021
% this function is for select the needed data segment
% Input:
%   x: a line vector dataset
%   t0: the start timepoints vector of segments you want to select.
%   Fs: sampling rate
%   dur: duration that you want to select from x (s)
% Output:
%   seg: the segmented data.

    dt = 1/Fs;
    seg = zeros(length(t0),dur*Fs); 
    for m = 1:size(seg,1)
        seg(m,:) = x(round((t0(m)+dt)*Fs):round((t0(m)+dur)*Fs));
    end
end