% this function is for calculate coherence by multi taper method
% Choose time-bandwidth product of 20.(10s*4hz/2), freq resolution is 4hz
function [C,phi,S12,S1,S2,f]=cohrmt(E1,E2,Fs,dur)

    T = dur;             %Define total duration of data,
    dt = 1/Fs;           %... sampling interval,
    N = T/dt;            %... and # pts in data.

    TW = T*4/2;				%Choose time-bandwidth product of 20.(10s*4hz/2), freq resolution is 4hz
    ntapers = 2*TW-1;		%...which sets the # of tapers.
    params.Fs = Fs;		%Define frequency resolution,
    params.tapers = [TW,ntapers];%...time-band product,# tapers.
    params.pad = -1;		%Specify no zero padding,
                            %... and compute the coherence.
    [C,phi,S12,S1,S2,f]=coherencyc(E1, E2, params);

end