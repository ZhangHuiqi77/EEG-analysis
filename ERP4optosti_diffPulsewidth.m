%% plot Event Related Potential (ERP)
% read data
hinfo = hdf5info('20200313_A79_CHETA_sti_M1R_rec_M1R1_20HZ_20ms_20mw_30s_180s_times3_0001');
dset = hdf5read(hinfo.GroupHierarchy.Groups(2).Datasets(1));
dset = dset';
dset_1 = double(dset(1,:)); % M1
dset_sti = double(dset(3,:)); % sti tag


% define sampling parameters
Fs = 10000;
gain = 10000;
dt = 1/Fs;
t = [0:dt:(size(dset_1,2)-1)*dt]; % creat t as x-axis
N = length(dset_1);
T = N*dt;
df = 1/max(T);
fNQ = 1/dt/2;


% define opto sti parameters
Fsti = 20; % stimuli frequency %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pwidth = 20; % stimuli pulsewidth %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pwidth = pwidth/1000;
dur = 30; % duration [s] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delay = 180 ;% delay between two stimulations [s]%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cyc = 3; % sti cycles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% locate the Time for the light on and creat stimuli events
[pks,locs] = findpeaks(dset_sti,'MinPeakHeight',14000,'MinPeakDistance',pwidth*Fs+10);
% ERP0 = zeros(length(locs),pwidth/1000*Fs*2);% prepare a zero matrix for fill in stimuli events
ERP0 = zeros(length(locs),20/1000*Fs*2);% prepare a zero matrix for fill in stimuli events, using max pulsewidth 20
index = int64(locs);
for k = 1:length(locs)
    ERP0(k,:) = dset_1(index(k):index(k)+20/1000*Fs*2-1);
end
ERP1202020 = ERP0/gain*1000; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate peak-to-peak amptitude
mnERP = mean(ERP1202020,1); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p2pA = max(mnERP(:)) - min(mnERP(:))
p2pA = peak2peak(mnERP)
[M I] = min(mnERP);
[M A] = max(mnERP);
latency = I*dt*1000 % [ms] latency = X_lowest-loc_on
% duration = A*dt*1000
a = find(abs(mnERP-(-100))<2);
duration = a*dt*1000
%% plot ERP
figure
t3 = [0:dt:20/1000*2-dt];
plot(t3, mean(ERP0200520,1),'LineWidth',2) % the first 0 means no opsin, then 05-5mw, 10-10mw, 20-20mw
hold on
plot(t3, mean(ERP0201020,1),'LineWidth',2)
plot(t3, mean(ERP0202020,1),'LineWidth',2)
plot(t3, mean(ERP1200520,1),'LineWidth',2)
plot(t3, mean(ERP1201020,1),'LineWidth',2)
plot(t3, mean(ERP1202020,1),'LineWidth',2)
xlabel('Time [s]')
ylabel('Voltage [mV]')
ylim([-500 200])
title('LFP trace under 473nm light stimulation in right M1- 20hz 5/10/20ms 20mw') 
legend('no opsin - 5ms','no opsin - 10ms','no opsin - 20ms','ChETA - 5ms','ChETA - 10ms','ChETA - 20ms','Location', 'SouthEast') %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot([0.01,0.01],[0 -500],'k','LineWidth',2)
legend boxoff
