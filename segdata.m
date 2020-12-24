% this function is for select the needed data segment
% x is the 1-line dataset
% t0 is the start timepoints of every segment you want to select, 
% it should be a 1-dim list
% Fs is the sampling rate, which always has been defined at the start of
% the analysis
% dur is the duration that you want to select from x (s)

function [seg]= segdata(x,t0,Fs,dur)
    dt = 1/Fs;
    seg = zeros(length(t0),dur*Fs); 
    for m = 1:size(seg,1)
        seg(m,:) = x(((t0(m)+dt)/dt):(t0(m)+dur)/dt);
    end
end