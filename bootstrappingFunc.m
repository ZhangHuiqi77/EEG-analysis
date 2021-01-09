function [mn,ciL,ciU] = bootstrappingFunc(x,nb)
% Huiqi
% 02/01/2021
% This function generates bootstraping data with confidential interval.
% Input
%   x: data, include many trials, each row represent a trial
%   nb: number of bootstrapping application, which will generate nb ERP0.
% Output
%   mn: ERP of data(no bootstrapping application)
%   ciL: low limitation of confidential interval for bootstrapped ERPs,
%   2.5% for unilateral
%   ciU: Up limitation of confidential interval for bootstrapped ERPs, 
%   97.5% for unilateral

    ntrials = size(x,1);
    ERP0 = zeros(3000,size(x,2));       %Create empty ERP variable.
    for k=1:nb                          %For each resampling,
        i=randsample(ntrials,ntrials,1);%... choose the trials,
        EEG0 = x(i,:);                  %... create resampled EEG,
        ERP0(k,:) = mean(EEG0,1);		%... save resampled ERP.
    end

    sERP0=sort(ERP0);                   %Sort each column of resampled ERP.
    ciL  =sERP0(0.025*size(ERP0,1),:);  %Determine lower CI.
    ciU  =sERP0(0.975*size(ERP0,1),:);  %Determine upper CI.
    mn = mean(x,1);                     %Determine ERP,
end