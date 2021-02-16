function xepoch = data4epoch(x,t0,Fs,dur,sec)
% huiqi
% 02/01/2021
% called function: segdata
% This function is for preparing data that is consistent of many epoches
% USAGE: xepoch = data4epoch(x,t0,Fs,dur,sec)
% Input:
%   x: data, should be a line vector
%   t0: start time to segment data, should be a vector
%   Fs: sampling frequency
%   dur: the length for segment data, should be the same length as t0.
%   sec: the length for each epoch.
% Output:
%   xepoch: the generated data
    N = sum(dur)/sec;            %Define number that all seg can provide
    xepoch = zeros(N,sec*Fs);    % reserve structure
    nsum = 0;                    %Define a inter var.
    
    for m = 1:length(t0)
        n = dur(m)/sec;          % Define number that single seg can provide
        seg = segdata(x,t0(m),Fs,dur(m));    % extract move data from x
        xepoch((nsum+1:nsum+n),:) = reshape(seg,[sec*Fs,n])'; % reshape the data to sec in an epoch
        nsum = nsum + n;
    end
end