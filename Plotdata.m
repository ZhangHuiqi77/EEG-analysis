%% HS3

% [head, dset] = edfread('20201204_GL8213_PTZsc80_3rdHS_APCR1_APCL2_dHCR3_S1L4.edf');
% dset_1 = dset(1,:); % Ch1
% dset_2 = dset(2,:); % Ch2
% dset_3 = dset(3,:); % Ch3
% dset_4 = dset(4,:); % Ch4
% dset_5 = dset(5,:); % volocity
%% HS2
[head, dset] = edfread('20210113_GL8661_KA3.5M_DREADD_day6_CNO_day_APCR1_dHCR2_S1R3_S1L4.edf');
dset_1 = touv(dset(1,:)); %
dset_2 = touv(dset(2,:)); %
dset_3 = touv(dset(3,:)); %
dset_4 = touv(dset(4,:)); %
dset_5 = dset(5,:);       % volocity, in raw data structure

ch1 = 'APCR';               %Define channel name
ch2 = 'LSR';
ch3 = 'dHCR';
ch4 = 'S1L';

% define sampling parameters
Fs = 1000; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gain = 1; 
dt = 1/Fs; 

%% define the prepared dataset
x = dset_2; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = (dt:dt:(size(x,2))*dt);

figure('color','w');
tiledlayout(2,1);
ax1 = nexttile;
plot(t,x,'k')
title(['raw EEG- ',ch3])

ax2 = nexttile;
plot(t,dset_5)
title('acceleration')
xlabel('Time [s]')


%% set the range
linkaxes([ax1,ax2],'x')
ax1.XLim = [18894 18914];
ax1.YLim = [-600 600];
ax2.YLim = [-150 150];

