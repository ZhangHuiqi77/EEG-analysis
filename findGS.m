% this function is for finding epileptic spikes according to 2 standards
% Condition 1)freq >2 hz, so that two ajacent locs_spk should not bigger than 500( when the sampling rate is 1000)
% Condition 2) duration > 10s, a)it shuld have at least 20 continuous spks; and b)end - start > 10s

% we set 100 gs matrix and value in 10s for every gs, duration of every GS is saved before crop the data to 10 sec 


% the input parameters: x: data; Fs: sampling frequency;
% functions need to call: calbs
% output parameters: locgs: the numerical index of gs; gs: the data value
% of gs; gs_dur: duration of every gs; gs_num: the number of all gs,
% interval: interval between every two gs. 

function [locgs,gs,gs_dur,gs_num,interval]=findGS(x, Fs, bs)
    % find epileptic spikes
    [pks,locs] = findpeaks(-x,'MinPeakProminence',2*bs);
    
    % build up zero matrix to store variables  
    gs = zeros(100,10*Fs); % build gs structure that is sure enough to save every gs, here we use 100
    locgs = zeros(100,10*Fs); % build gs structure that is sure enough to save every gs location, here we use 100, the same as gs.
    locspk = zeros(1000,100000); % an inter variable, build a matrix to save location of spikes, the matrix is big enough, here we use length(x)
    gs_dur = zeros(100,1);
    
    % condition 1
    m = 1; 
    n = 1;
    for h = 1: (length(locs)-1) %split the spks to GSs
        if locs(h+1) - locs(h) < 500
            locspk(m,n) = locs(h+1);                                                                                                                                                                                                                                                                                                                                 
            n = n+1;
        else
            m = m+1;
            n = 1;
            locspk(m,n) = locs(h+1);
            n = n+1;
        end
    end
    
    % condition 2
    h = 1;  % initiate the number of GS. 
    for k = 1:size(locspk,1)
        locspk_non0 = nonzeros(locspk(k,:));
        if length(locspk_non0) > 20
            dur = (locspk_non0(end) - locspk_non0(1))/Fs;
            if dur > 10
                locgs(h,1:length(locspk_non0)) = locspk_non0';
                gs(h,1:10*Fs) = x(locspk_non0:(locspk_non0+10*Fs-1)); % assign 10s values to gs
                gs_dur(h,1) = dur; % save the duration of every GS in gs_duration, although only the 10s are saved in the variable "gs"
                h = h+1;
            end
        end
    end
    
    % the number of gs
    gs_num = length(nonzeros(gs(:,1)));
    
    % calculate interval
    interval = zeros(size(locgs,1)-1,1);
    h = 1;
    for h = 1: (size(locgs,1)-1)
        interval(h,1) = (locgs(h+1,1)-max(locgs(h,:)))/1000;
    end
end



% check variables
% gs -- generalized seizure voltage
% locs2 -- gs index
% locs -- the number of spikes
% gs_dur -- duration of every GS
% interval -- interval between every two GS.
