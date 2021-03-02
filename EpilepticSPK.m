% find epileptic spikes
%read data- edf
[head, dset] = edfread('20201204_GL8213_PTZsc80_3rdHS_APCR1_APCL2_dHCR3_S1L4.edf');
dset_1 = touv(dset(1,:)); 
dset_2 = touv(dset(2,:)); 
dset_3 = touv(dset(3,:)); 
dset_4 = touv(dset(4,:)); 
dset_5 = dset(5,:);
% define sampling parameters
Fs = 1000; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gain = 1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = 1/Fs; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define the prepared dataset
x = dset_3; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = (dt:dt:(size(x,2))*dt);

figure('color','w')
% hold on
plot(t,x,'k')
% plot(t,dset_2-1000,'k')
% plot(t,dset_3-1000*2,'k')
% plot(t,dset_4-1000*3,'k')
% ylim([-10000,1000])
title('raw EEG-APCR') %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plot all 4 channels
xlimit = [1500 1540]; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ylimit = [-2000 2000]; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(5,1,1)
plot(t,dset_1,'k')
% xlabel('Time [s]')
ylabel('Voltage [\mu V]')
% xlim([xlimit(1) xlimit(2)]);
ylim([ylimit(1) ylimit(2)]);
title('raw EEG - APCR')

subplot(5,1,2)
plot(t,dset_2,'k')
% xlabel('Time [s]')
ylabel('Voltage [\mu V]')
% xlim([xlimit(1) xlimit(2)]);
ylim([ylimit(1) ylimit(2)]);
title('raw EEG - APCL')

subplot(5,1,3)
plot(t,dset_3,'k')
% xlabel('Time [s]')
ylabel('Voltage [\mu V]')
% xlim([xlimit(1) xlimit(2)]);
ylim([ylimit(1) ylimit(2)]);
title('raw EEG - S1R')

subplot(5,1,4)
plot(t,dset_4,'k')
% xlabel('Time [s]')
ylabel('Voltage [\mu V]')
% xlim([xlimit(1) xlimit(2)]);
ylim([ylimit(1) ylimit(2)]);
title('raw EEG - S1L')

subplot(5,1,5)
plot(t,dset_5)
xlabel('Time [s]')
% xlim([xlimit(1) xlimit(2)]);
ylabel('acceleration')
% define baseline
% define pre-ictal to calculate baseline
% baseline is an interval that 98% data are dropped in.

t0 = (0);
x_bs = segdata(x,t0,Fs,600);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%0-600s as bs stage
bs = calbs(x_bs);

% plot epileptic peaks(p2pA > 3*bs_p2p)
[pks,locs] = findpeaks(-x,'MinPeakProminence',2*bs);
figure
plot(t,-x,'k',t(locs),pks,'or')

% find gs 1) Freq > 2hz, 2) dur > 10s
[locgs,gs,gs_dur,gs_num,interval]=findGS(x, Fs,bs);

% define t_pregs, t_gs, 10s       for 1 animal
t_ini = locgs(1:gs_num,1)/Fs;
x_pregs = segdata(x,t_ini-10,Fs,10);
x_gs = segdata(x,t_ini,Fs,10);

% concatenates data in the same group

x_pregs_ctrl1 = [x_pregs_ctrl1;x_pregs];
x_gs_ctrl1 = [x_gs_ctrl1;x_gs];

