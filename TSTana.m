% This script is for TST analysis.
% this script provide common processing methods for LFP data
% Called hand-made functions: touv, segdata,data4epoch, shuffleFunc,
% bootstrappingFunc, trialAveSpectrum, spectrogramFunc,cohrFunc

%% plot the data
% read data- edf
[head, dset] = edfread('20201103_GL8214_TST_box2_APCR1_LSR2_dHCR3_S1L4.edf');
dset_1 = touv(dset(1,:));   % APCR %%%%%%%%%%%%%%%%%%%%%%%%%
dset_2 = touv(dset(2,:));   % APCL %%%%%%%%%%%%%%%%%%%%%%%%
dset_3 = touv(dset(3,:));   % dHCR %%%%%%%%%%%%%%%%%%%%%%%%%
dset_4 = touv(dset(4,:));   % S1L %%%%%%%%%%%%%%%%%%%%%%%%%
dset_5 = dset(5,:);
dset = [dset_1;dset_2;dset_3;dset_4;dset_5]; % rebuild dset after touv transform.

ch1 = 'APCR';               %Define channel name
ch2 = 'LSR';
ch3 = 'dHCR';
ch4 = 'S1L';
% define sampling parameters
Fs = 1000; 
dt = 1/Fs; 
t = (dt:dt:(size(dset_1,2))*dt);

figure;plot(t,dset_5)
x4bs = segdata(dset_5,193,1000,10);  %Define data for baseline cal by eye.


% plot all 4 channels
figure
tiledlayout(5,1);
ax1 = nexttile;
plot(t,dset_1,'k')
% title('raw LFP - APCR')

ax2 = nexttile;
plot(t,dset_2,'k')
% title('raw LFP - APCL')

ax3 = nexttile;
plot(t,dset_3,'k')
% title('raw LFP - dHCR')

ax4 = nexttile;
plot(t,dset_4,'k')
% title('raw EEG - S1R')
% xlabel('Time [s]')
ylabel('Voltage [\mu V]')

ax5 = nexttile;
plot(t,dset_5)
% title('acceleration')
xlabel('Time [s]')

linkaxes([ax1,ax2,ax3,ax4,ax5],'x')
ax1.XLim = [0 30]; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
linkaxes([ax1,ax2,ax3,ax4],'y')
ax1.YLim = [-300 300];  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% preliminary spectrum 
% this is a roughly check of the total timeline. 
x = dset_3; %%%%%%%%%%%%%%%%%%%%%%%%% Define interested Ch
t = (dt:dt:(size(x,2))*dt);
N = length(x);
T = t(end);

x = hann(N).*transpose(x);      %Multiply data by Hanning window.
xf = fft(x-mean(x));			%Compute Fourier transform of x.
Sxx = 2*dt^2/T *(xf.*conj(xf));	%Compute the power.
Sxx = Sxx(1:N/2+1);				%Ignore negative frequencies.

df = 1/max(T);					%Define frequency resolution.
fNQ = 1/dt/2;					%Define Nyquist frequency.
faxis = (0:df:fNQ);				%Construct frequency axis.
figure
plot(faxis, 10*log10(Sxx))		%Plot power vs frequency.
xlim([0 100]); ylim([-60 60])	%Set frequency & power ranges.
xlabel('Frequency [Hz]')		%Label axes.
ylabel('Power [ mV^2/Hz]')

% figure;plot(faxis,Sxx)
% xlim([0 200])
% ylim([0 100])
 
%% split to bands and link axis
x = dset_3; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = t(2)-t(1);				%Define the sampling interval.
Fs = 1/dt;					%Define the sampling frequency.
fNQ = Fs/2;					%Define the Nyquist frequency.
							%For low frequency interval,
Wn = [5,12]/fNQ;				%...set the passband, Theta
n  = 100;					%...and filter order,
b  = fir1(n,Wn);			%...build bandpass filter.
Vlo = filtfilt(b,1,x);      %...and apply filter.
							%For high frequency interval,
Wn = [30,50]/fNQ;			%...set the passband, low Gamma
n  = 100;					%...and filter order,
b  = fir1(n,Wn);			%...build bandpass filter.
Vhi = filtfilt(b,1,x);  	%...and apply filter.


figure;tiledlayout(4,1);
ax1 = nexttile;
plot(t,x,'k')

ax2 = nexttile;
plot(t,Vlo,'k')

ax3 = nexttile;
plot(t,Vhi,'k')

ax4 = nexttile;
plot(t,dset_5)

linkaxes([ax1,ax2,ax3,ax4],'x')
ax1.XLim = [195 200]; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linkaxes([ax1,ax2],'y')
ax1.YLim = [-200 200];

%% extract the data automatically
% our aim is to extract epoches of move and unmove stage, 
% ...according to channel 5, movement.
x = dset_5-mean(dset_5);        %Define the movement channel.
Fs = 1000;                      %Define sampling frequency.
dt = 1/Fs;                      %Define time interval for sampling.
t = (dt:dt:(size(x,2))*dt);     %Define time variable. 

bs = calbs(x4bs);               %Calculate basiline- p2p difference.
[pks_m,locs_m] = findpeaks(x,t,'MinPeakProminence',3*bs);   
                                %...Define movement threshold, 3*bs.

figure;
plot(t,x)                   %...plot movement.
hold on
plot(locs_m,pks_m,'o')      %...plot the found move peaks
xlabel('time [s]')
hold off

% detect the unmove time(> 1s)
difflocs_m = diff(locs_m);  %...calculate time interval of the move point.
figure;
plot(difflocs_m)
ylabel('moving time interval[s]')

thr = 1;    %Define a threshold for the unmove duration.
[pks,locs] = findpeaks(difflocs_m,'MinPeakHeight',thr); 
hold on
plot(locs,pks,'o')
t0um = locs_m(locs);    %start time loc of unmove.
t0m = locs_m(locs+1);   %start time loc of move.
t0m = [0 t0m];          %add 0 to starting point of m.

figure;                 %plot and have a check
plot(t,x)
hold on
for n = 1:length(t0m)
    plot([t0m(n) t0m(n)],[-200 200],'r');
    hold on
end

for n = 1:length(t0um)
    plot([t0um(n) t0um(n)],[-200 200],'k');
    hold on
end

%% prepare data: move epoch
t0mm = t0m(1:(length(t0m)-1)); %Define t0m as m-um loop. so rm the last one.
t0mu = t0m(2:end);             %Define t0m as um-m loop.
mdur = t0um-t0mm;              % calculate the duration of every movement. 
mdur = floor(mdur);            % round the mdur to integers.
ind = find(mdur==0);           % find the index of mdur is 0, 
% ... which is too short and won't be selected in the later analysis.
t0mm(ind) = [];             % remove t0mm that dur is 0.
t0mmu = t0um;               % prepare a new t0um, corresponding to t0mm
t0mmu(ind)=[];              % remove coresponding t0um.
mdur = mdur(mdur>0);        % remain the mdur >=1, and it is integer.

figure;plot(t,x)    % ...plot and eye checking time for move start.
hold on
for n = 1:length(t0mm)
    plot([t0mm(n) t0mm(n)],[-200 200],'r');
    hold on
end
for n = 1:length(t0mmu)
    plot([t0mmu(n) t0mmu(n)],[-200 200],'k');
    hold on
end
hold off

moveEpoch =  data4epoch(x,t0mm,1000,mdur,1);    
moveEpoch1 = data4epoch(dset_1,t0mm,1000,mdur,1);   % xepoch = data4epoch(x,t0,Fs,dur,sec),Ch1
moveEpoch2 = data4epoch(dset_2,t0mm,1000,mdur,1);   % xepoch = data4epoch(x,t0,Fs,dur,sec),Ch2
moveEpoch3 = data4epoch(dset_3,t0mm,1000,mdur,1);   % xepoch = data4epoch(x,t0,Fs,dur,sec),Ch3
moveEpoch4 = data4epoch(dset_4,t0mm,1000,mdur,1);   % xepoch = data4epoch(x,t0,Fs,dur,sec),Ch4


%% prepare data: unmove epoch
umdur = t0mu-t0um;              % calculate the duration of every movement. 
umdur = floor(umdur);            % round the mdur to integers.

figure;plot(t,dset_5)    % ...plot and eye checking time for move start.
hold on
for n = 1:length(t0um)
    plot([t0um(n) t0um(n)],[-200 200],'k');
    hold on
end
for n = 1:length(t0mu)
    plot([t0mu(n) t0mu(n)],[-200 200],'r');
    hold on
end
hold off

unmoveEpoch =  data4epoch(x,t0um,1000,umdur,1);
unmoveEpoch1 = data4epoch(dset_1,t0um,1000,umdur,1);   % xepoch = data4epoch(x,t0,Fs,dur,sec),Ch1
unmoveEpoch2 = data4epoch(dset_2,t0um,1000,umdur,1);   % xepoch = data4epoch(x,t0,Fs,dur,sec),Ch2
unmoveEpoch3 = data4epoch(dset_3,t0um,1000,umdur,1);   % xepoch = data4epoch(x,t0,Fs,dur,sec),Ch3
unmoveEpoch4 = data4epoch(dset_4,t0um,1000,umdur,1);   % xepoch = data4epoch(x,t0,Fs,dur,sec),Ch4

%% prepare data ( >5s m2um um2m)
t0mm = t0m(1:(length(t0m)-1)); %Define t0m as m-um loop. so rm the last one.
t0mu = t0m(2:end);             %Define t0m as um-m loop.
mdur = t0um-t0mm;              % calculate the duration of every movement.
umdur = t0mu-t0um;              % calculate the duration of every unmovement.

nloop = length(umdur);     % m-um loop
m2um = zeros(nloop,10*Fs);
m2um1 = zeros(nloop,10*Fs);
m2um2 = zeros(nloop,10*Fs);
m2um3 = zeros(nloop,10*Fs);
m2um4 = zeros(nloop,10*Fs);

um2m = zeros(nloop,10*Fs);
um2m1 = zeros(nloop,10*Fs);
um2m2 = zeros(nloop,10*Fs);
um2m3 = zeros(nloop,10*Fs);
um2m4 = zeros(nloop,10*Fs);

for m = 1:nloop
   if mdur(m)>5 && umdur(m)>5       % m-um loop both > 5s
       m2um(m,:) = segdata(dset_5,t0um(m)-5,Fs,10);      % m2um ch5
       m2um1(m,:) = segdata(dset_1,t0um(m)-5,Fs,10);% m2um ch1
       m2um2(m,:) = segdata(dset_2,t0um(m)-5,Fs,10);% m2um ch2
       m2um3(m,:) = segdata(dset_3,t0um(m)-5,Fs,10);% m2um ch3
       m2um4(m,:) = segdata(dset_4,t0um(m)-5,Fs,10);% m2um ch4
   end
end
for m = 1:nloop-1
    if umdur(m)>5 && mdur(m+1)>5    % um-m loop both > 5s
        um2m(m,:) = segdata(dset_5,t0mu(m)-5,Fs,10);           % um2m ch5
        um2m1(m,:) = segdata(dset_1,t0mu(m)-5,Fs,10);     % um2m ch1
        um2m2(m,:) = segdata(dset_2,t0mu(m)-5,Fs,10);     % um2m ch2
        um2m3(m,:) = segdata(dset_3,t0mu(m)-5,Fs,10);     % um2m ch3
        um2m4(m,:) = segdata(dset_4,t0mu(m)-5,Fs,10);     % um2m ch4
    end
end
m2um = m2um(any(m2um,2),any(m2um,1));       % remain nonzero lines.
m2um1 = m2um1(any(m2um1,2),any(m2um1,1));   % remain nonzero lines.
m2um2 = m2um2(any(m2um2,2),any(m2um2,1));   % remain nonzero lines.
m2um3 = m2um3(any(m2um3,2),any(m2um3,1));   % remain nonzero lines.
m2um4 = m2um4(any(m2um4,2),any(m2um4,1));   % remain nonzero lines.

um2m = um2m(any(um2m,2),any(um2m,1));       % remain nonzero lines.
um2m1 = um2m1(any(um2m1,2),any(um2m1,1));   % remain nonzero lines.
um2m2 = um2m2(any(um2m2,2),any(um2m2,1));   % remain nonzero lines.
um2m3 = um2m3(any(um2m3,2),any(um2m3,1));   % remain nonzero lines.
um2m4 = um2m4(any(um2m4,2),any(um2m4,1));   % remain nonzero lines.
%% 1: ERP 
x = m2um1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = (dt:dt:10);
figure;
subplot(4,2,1);
hold on
for m = 1:size(x,1)
    plot(t,x,'r')
end
hold on
plot(t,mean(x),'LineWidth',2,'color','k')
plot([5 5],[-200 200],'LineWidth',1,'color','w');
plot([0 10],[0 0],'LineWidth',1,'color','w')
ylim([-50 50])
title(['m2um-',ch1])
hold off

x = um2m1; %%%%%%%%%%%%%%%%%%%%%%%%%%%
t = (dt:dt:10);
subplot(4,2,2)
hold on
for m = 1:size(x,1)
    plot(t,x,'r')
end
plot(t,mean(x),'LineWidth',2,'color','k')
plot([5 5],[-150 200],'LineWidth',1,'color','w');
plot([0 10],[0 0],'LineWidth',1,'color','w')
ylim([-50 50])
title(['um2m-',ch1])%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold off

x = m2um2; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = (dt:dt:10);
subplot(4,2,3);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
for m = 1:size(x,1)
    plot(t,x,'r')
end
hold on
plot(t,mean(x),'LineWidth',2,'color','k')
plot([5 5],[-200 200],'LineWidth',1,'color','w');
plot([0 10],[0 0],'LineWidth',1,'color','w')
ylim([-50 50])
title(['m2um-',ch2])%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold off

x = um2m2; %%%%%%%%%%%%%%%%%%%%%%%%%%%
t = (dt:dt:10);
subplot(4,2,4)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
for m = 1:size(x,1)
    plot(t,x,'r')
end
plot(t,mean(x),'LineWidth',2,'color','k')
plot([5 5],[-150 200],'LineWidth',1,'color','w');
plot([0 10],[0 0],'LineWidth',1,'color','w')
ylim([-50 50])
title(['um2m-',ch2])%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold off

x = m2um3; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = (dt:dt:10);
subplot(4,2,5);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
for m = 1:size(x,1)
    plot(t,x,'r')
end
hold on
plot(t,mean(x),'LineWidth',2,'color','k')
plot([5 5],[-200 200],'LineWidth',1,'color','w');
plot([0 10],[0 0],'LineWidth',1,'color','w')
ylim([-50 50])
title(['m2um-',ch3])%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold off

x = um2m3; %%%%%%%%%%%%%%%%%%%%%%%%%%%
t = (dt:dt:10);
subplot(4,2,6)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
for m = 1:size(x,1)
    plot(t,x,'r')
end
plot(t,mean(x),'LineWidth',2,'color','k')
plot([5 5],[-150 200],'LineWidth',1,'color','w');
plot([0 10],[0 0],'LineWidth',1,'color','w')
ylim([-50 50])
title(['um2m-',ch3])%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold off

x = m2um4; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = (dt:dt:10);
subplot(4,2,7);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
for m = 1:size(x,1)
    plot(t,x,'r')
end
hold on
plot(t,mean(x),'LineWidth',2,'color','k')
plot([5 5],[-200 200],'LineWidth',1,'color','w');
plot([0 10],[0 0],'LineWidth',1,'color','w')
ylim([-50 50])
title(['m2um-',ch4])%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold off

x = um2m4; %%%%%%%%%%%%%%%%%%%%%%%%%%%
t = (dt:dt:10);
subplot(4,2,8)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
for m = 1:size(x,1)
    plot(t,x,'r')
end
plot(t,mean(x),'LineWidth',2,'color','k')
plot([5 5],[-150 200],'LineWidth',1,'color','w');
plot([0 10],[0 0],'LineWidth',1,'color','w')
ylim([-50 50])
title(['um2m-',ch4])%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold off
% the result of oringinal selected trials(4 and 5 trials) shows a little mess and 
% seems no significance in the m2um data, some significant change in um2m
% data. Anyway, bootstrapping give a more convinced result with more repeat

%% 2: plot bootstrapping ERP of m2um and um2m
dt = 0.001;
t = (-5:dt:5-dt);
% Ch1
figure;
subplot(4,2,1)
[mn,ciL,ciU] = bootstrappingFunc(m2um1,3000);  % apply bootstrapping for 3000 times
plot(t,ciL,'b')                             %... and plot lower CI,
hold on
plot(t,ciU,'b')                             %... and upper CI.
plot(t, mn, 'LineWidth', 3,'color','k')     %... and plot it.
ylabel('Voltage [\mu V]')       %Label the y-axis as voltage.
xlabel('Time [s]')              %Label the x-axis as time.
title(['ERP of m2um-',ch1, ' with bootstrapping CI'])
plot([0 0],[-80 80],'LineWidth',1,'color','w');
plot([-5 5],[0 0],'LineWidth',1,'color','w')
hold off

subplot(4,2,2)
[mn,ciL,ciU] = bootstrappingFunc(um2m1,3000);  % apply bootstrapping for 3000 times
plot(t,ciL,'b')                     %... and plot lower CI,
hold on
plot(t,ciU,'b')                     %... and upper CI.
plot(t, mn, 'LineWidth', 3,'color','k')    %... and plot it.
ylabel('Voltage [\mu V]')       %Label the y-axis as voltage.
xlabel('Time [s]')              %Label the x-axis as time.
title(['ERP of um2m-',ch1,' with bootstrapping CI'])
plot([0 0],[-80 80],'LineWidth',1,'color','w');
plot([-5 5],[0 0],'LineWidth',1,'color','w')
hold off

% Ch2
subplot(4,2,3)
[mn,ciL,ciU] = bootstrappingFunc(m2um2,3000);  % apply bootstrapping for 3000 times
plot(t,ciL,'b')                     %... and plot lower CI,
hold on
plot(t,ciU,'b')                     %... and upper CI.
plot(t, mn, 'LineWidth', 3,'color','k')    %... and plot it.
ylabel('Voltage [\mu V]')       %Label the y-axis as voltage.
xlabel('Time [s]')              %Label the x-axis as time.
title(['ERP of m2um-',ch2, ' with bootstrapping CI'])
plot([0 0],[-80 80],'LineWidth',1,'color','w');
plot([-5 5],[0 0],'LineWidth',1,'color','w')
hold off

subplot(4,2,4)
[mn,ciL,ciU] = bootstrappingFunc(um2m2,3000);  % apply bootstrapping for 3000 times
plot(t,ciL,'b')                     %... and plot lower CI,
hold on
plot(t,ciU,'b')                     %... and upper CI.
plot(t, mn, 'LineWidth', 3,'color','k')    %... and plot it.
ylabel('Voltage [\mu V]')       %Label the y-axis as voltage.
xlabel('Time [s]')              %Label the x-axis as time.
title(['ERP of um2m-',ch2,' with bootstrapping CI'])
plot([0 0],[-80 80],'LineWidth',1,'color','w');
plot([-5 5],[0 0],'LineWidth',1,'color','w')
hold off

% Ch3
subplot(4,2,5)
[mn,ciL,ciU] = bootstrappingFunc(m2um3,3000);  % apply bootstrapping for 3000 times
plot(t,ciL,'b')                     %... and plot lower CI,
hold on
plot(t,ciU,'b')                     %... and upper CI.
plot(t, mn, 'LineWidth', 3,'color','k')    %... and plot it.
ylabel('Voltage [\mu V]')       %Label the y-axis as voltage.
xlabel('Time [s]')              %Label the x-axis as time.
title(['ERP of m2um-',ch3, ' with bootstrapping CI'])
plot([0 0],[-80 80],'LineWidth',1,'color','w');
plot([-5 5],[0 0],'LineWidth',1,'color','w')
hold off

subplot(4,2,6)
[mn,ciL,ciU] = bootstrappingFunc(um2m3,3000);  % apply bootstrapping for 3000 times
plot(t,ciL,'b')                     %... and plot lower CI,
hold on
plot(t,ciU,'b')                     %... and upper CI.
plot(t, mn, 'LineWidth', 3,'color','k')    %... and plot it.
ylabel('Voltage [\mu V]')       %Label the y-axis as voltage.
xlabel('Time [s]')              %Label the x-axis as time.
title(['ERP of um2m-',ch3,' with bootstrapping CI'])
plot([0 0],[-80 80],'LineWidth',1,'color','w');
plot([-5 5],[0 0],'LineWidth',1,'color','w')
hold off

% Ch4
subplot(4,2,7)
[mn,ciL,ciU] = bootstrappingFunc(m2um4,3000);  % apply bootstrapping for 3000 times
plot(t,ciL,'b')                     %... and plot lower CI,
hold on
plot(t,ciU,'b')                     %... and upper CI.
plot(t, mn, 'LineWidth', 3,'color','k')    %... and plot it.
ylabel('Voltage [\mu V]')       %Label the y-axis as voltage.
xlabel('Time [s]')              %Label the x-axis as time.
title(['ERP of m2um-',ch4, ' with bootstrapping CI'])
plot([0 0],[-80 80],'LineWidth',1,'color','w');
plot([-5 5],[0 0],'LineWidth',1,'color','w')
hold off

subplot(4,2,8)
[mn,ciL,ciU] = bootstrappingFunc(um2m4,3000);  % apply bootstrapping for 3000 times
plot(t,ciL,'b')                     %... and plot lower CI,
hold on
plot(t,ciU,'b')                     %... and upper CI.
plot(t, mn, 'LineWidth', 3,'color','k')    %... and plot it.
ylabel('Voltage [\mu V]')       %Label the y-axis as voltage.
xlabel('Time [s]')              %Label the x-axis as time.
title(['ERP of um2m-',ch4,' with bootstrapping CI'])
plot([0 0],[-80 80],'LineWidth',1,'color','w');
plot([-5 5],[0 0],'LineWidth',1,'color','w')
hold off

%% 3: plot bootstrapping ERP of m2um & um2m with shuffled signals
dt = 0.001;
t = (-5:dt:5-dt);
% Ch1
% prepare shuffled data
xs = shuffleFunc(m2um1,1);                  % generate shuffled data
[mn,ciL,ciU] = bootstrappingFunc(xs,3000);  % apply bootstrapping to shuffled data for 3000 times
figure;
subplot(4,2,1)
plot(t,ciL,'color',[0.8 0.8 0.8])           %... and plot lower CI,
hold on
plot(t,ciU,'color',[0.8 0.8 0.8])           %... and upper CI.
plot(t, mn, 'LineWidth', 3,'color', [0.5 0.5 0.5])    %... and plot it.
% bootstrapping for target data
[mn,ciL,ciU] = bootstrappingFunc(m2um1,3000);  % apply bootstrapping for 3000 times
plot(t,ciL,'b')                     %... and plot lower CI,
hold on
plot(t,ciU,'b')                     %... and upper CI.
plot(t, mn, 'LineWidth', 3,'color','k')    %... and plot it.
ylabel('Voltage [\mu V]')       %Label the y-axis as voltage.
xlabel('Time [s]')              %Label the x-axis as time.
title(['ERP of m2um-',ch1,' with bootstrapping CI and shuffled signals'])
plot([0 0],[-80 80],'LineWidth',1,'color','w');
plot([-5 5],[0 0],'LineWidth',1,'color','w')
hold off

% bootstrapping for shuffled data
xs = shuffleFunc(um2m1,1);                  % generate shuffled data
[mn,ciL,ciU] = bootstrappingFunc(xs,3000);  % apply bootstrapping to shuffled data for 3000 times
subplot(4,2,2)
plot(t,ciL,'color',[0.8 0.8 0.8])           %... and plot lower CI,
hold on
plot(t,ciU,'color',[0.8 0.8 0.8])           %... and upper CI.
plot(t, mn, 'LineWidth', 3,'color', [0.5 0.5 0.5])    %... and plot it.
% bootstrapping for target data
[mn,ciL,ciU] = bootstrappingFunc(um2m1,3000);  % apply bootstrapping for 3000 times
plot(t,ciL,'b')                     %... and plot lower CI,
hold on
plot(t,ciU,'b')                     %... and upper CI.
plot(t, mn, 'LineWidth', 3,'color','k')    %... and plot it.
ylabel('Voltage [\mu V]')       %Label the y-axis as voltage.
xlabel('Time [s]')              %Label the x-axis as time.
title(['ERP of um2m-',ch1,' with bootstrapping CI and shuffled signals'])
plot([0 0],[-80 80],'LineWidth',1,'color','w');
plot([-5 5],[0 0],'LineWidth',1,'color','w')
hold off

% Ch2
% prepare shuffled data
xs = shuffleFunc(m2um2,1);                  % generate shuffled data
[mn,ciL,ciU] = bootstrappingFunc(xs,3000);  % apply bootstrapping to shuffled data for 3000 times
subplot(4,2,3)
plot(t,ciL,'color',[0.8 0.8 0.8])           %... and plot lower CI,
hold on
plot(t,ciU,'color',[0.8 0.8 0.8])           %... and upper CI.
plot(t, mn, 'LineWidth', 3,'color', [0.5 0.5 0.5])    %... and plot it.
% bootstrapping for target data
[mn,ciL,ciU] = bootstrappingFunc(m2um2,3000);  % apply bootstrapping for 3000 times
plot(t,ciL,'b')                     %... and plot lower CI,
hold on
plot(t,ciU,'b')                     %... and upper CI.
plot(t, mn, 'LineWidth', 3,'color','k')    %... and plot it.
ylabel('Voltage [\mu V]')       %Label the y-axis as voltage.
xlabel('Time [s]')              %Label the x-axis as time.
title(['ERP of m2um-',ch2,' with bootstrapping CI and shuffled signals'])
plot([0 0],[-80 80],'LineWidth',1,'color','w');
plot([-5 5],[0 0],'LineWidth',1,'color','w')
hold off

% bootstrapping for shuffled data
xs = shuffleFunc(um2m2,1);                  % generate shuffled data
[mn,ciL,ciU] = bootstrappingFunc(xs,3000);  % apply bootstrapping to shuffled data for 3000 times
subplot(4,2,4)
plot(t,ciL,'color',[0.8 0.8 0.8])           %... and plot lower CI,
hold on
plot(t,ciU,'color',[0.8 0.8 0.8])           %... and upper CI.
plot(t, mn, 'LineWidth', 3,'color', [0.5 0.5 0.5])    %... and plot it.
% bootstrapping for target data
[mn,ciL,ciU] = bootstrappingFunc(um2m2,3000);  % apply bootstrapping for 3000 times
plot(t,ciL,'b')                     %... and plot lower CI,
hold on
plot(t,ciU,'b')                     %... and upper CI.
plot(t, mn, 'LineWidth', 3,'color','k')    %... and plot it.
ylabel('Voltage [\mu V]')       %Label the y-axis as voltage.
xlabel('Time [s]')              %Label the x-axis as time.
title(['ERP of um2m-',ch2,' with bootstrapping CI and shuffled signals'])
plot([0 0],[-80 80],'LineWidth',1,'color','w');
plot([-5 5],[0 0],'LineWidth',1,'color','w')
hold off

% Ch3
% prepare shuffled data
xs = shuffleFunc(m2um3,1);                  % generate shuffled data
[mn,ciL,ciU] = bootstrappingFunc(xs,3000);  % apply bootstrapping to shuffled data for 3000 times
subplot(4,2,5)
plot(t,ciL,'color',[0.8 0.8 0.8])           %... and plot lower CI,
hold on
plot(t,ciU,'color',[0.8 0.8 0.8])           %... and upper CI.
plot(t, mn, 'LineWidth', 3,'color', [0.5 0.5 0.5])    %... and plot it.
% bootstrapping for target data
[mn,ciL,ciU] = bootstrappingFunc(m2um3,3000);  % apply bootstrapping for 3000 times
plot(t,ciL,'b')                     %... and plot lower CI,
hold on
plot(t,ciU,'b')                     %... and upper CI.
plot(t, mn, 'LineWidth', 3,'color','k')    %... and plot it.
ylabel('Voltage [\mu V]')       %Label the y-axis as voltage.
xlabel('Time [s]')              %Label the x-axis as time.
title(['ERP of m2um-',ch3,' with bootstrapping CI and shuffled signals'])
plot([0 0],[-80 80],'LineWidth',1,'color','w');
plot([-5 5],[0 0],'LineWidth',1,'color','w')
hold off

% bootstrapping for shuffled data
xs = shuffleFunc(um2m3,1);                  % generate shuffled data
[mn,ciL,ciU] = bootstrappingFunc(xs,3000);  % apply bootstrapping to shuffled data for 3000 times
subplot(4,2,6)
plot(t,ciL,'color',[0.8 0.8 0.8])           %... and plot lower CI,
hold on
plot(t,ciU,'color',[0.8 0.8 0.8])           %... and upper CI.
plot(t, mn, 'LineWidth', 3,'color', [0.5 0.5 0.5])    %... and plot it.
% bootstrapping for target data
[mn,ciL,ciU] = bootstrappingFunc(um2m3,3000);  % apply bootstrapping for 3000 times
plot(t,ciL,'b')                     %... and plot lower CI,
hold on
plot(t,ciU,'b')                     %... and upper CI.
plot(t, mn, 'LineWidth', 3,'color','k')    %... and plot it.
ylabel('Voltage [\mu V]')       %Label the y-axis as voltage.
xlabel('Time [s]')              %Label the x-axis as time.
title(['ERP of um2m-',ch3,' with bootstrapping CI and shuffled signals'])
plot([0 0],[-80 80],'LineWidth',1,'color','w');
plot([-5 5],[0 0],'LineWidth',1,'color','w')
hold off

% Ch4
% prepare shuffled data
xs = shuffleFunc(m2um4,1);                  % generate shuffled data
[mn,ciL,ciU] = bootstrappingFunc(xs,3000);  % apply bootstrapping to shuffled data for 3000 times
subplot(4,2,7)
plot(t,ciL,'color',[0.8 0.8 0.8])           %... and plot lower CI,
hold on
plot(t,ciU,'color',[0.8 0.8 0.8])           %... and upper CI.
plot(t, mn, 'LineWidth', 3,'color', [0.5 0.5 0.5])    %... and plot it.
% bootstrapping for target data
[mn,ciL,ciU] = bootstrappingFunc(m2um4,3000);  % apply bootstrapping for 3000 times
plot(t,ciL,'b')                     %... and plot lower CI,
hold on
plot(t,ciU,'b')                     %... and upper CI.
plot(t, mn, 'LineWidth', 3,'color','k')    %... and plot it.
ylabel('Voltage [\mu V]')       %Label the y-axis as voltage.
xlabel('Time [s]')              %Label the x-axis as time.
title(['ERP of m2um-',ch4,' with bootstrapping CI and shuffled signals'])
plot([0 0],[-80 80],'LineWidth',1,'color','w');
plot([-5 5],[0 0],'LineWidth',1,'color','w')
hold off

% bootstrapping for shuffled data
xs = shuffleFunc(um2m4,1);                  % generate shuffled data
[mn,ciL,ciU] = bootstrappingFunc(xs,3000);  % apply bootstrapping to shuffled data for 3000 times
subplot(4,2,8)
plot(t,ciL,'color',[0.8 0.8 0.8])           %... and plot lower CI,
hold on
plot(t,ciU,'color',[0.8 0.8 0.8])           %... and upper CI.
plot(t, mn, 'LineWidth', 3,'color', [0.5 0.5 0.5])    %... and plot it.
% bootstrapping for target data
[mn,ciL,ciU] = bootstrappingFunc(um2m4,3000);  % apply bootstrapping for 3000 times
plot(t,ciL,'b')                     %... and plot lower CI,
hold on
plot(t,ciU,'b')                     %... and upper CI.
plot(t, mn, 'LineWidth', 3,'color','k')    %... and plot it.
ylabel('Voltage [\mu V]')       %Label the y-axis as voltage.
xlabel('Time [s]')              %Label the x-axis as time.
title(['ERP of um2m-',ch4,' with bootstrapping CI and shuffled signals'])
plot([0 0],[-80 80],'LineWidth',1,'color','w');
plot([-5 5],[0 0],'LineWidth',1,'color','w')
hold off

% shuffled data has no rhythm and almost cover the y-lim.

% the mean of dset_1 2 3 4 5 are in 0.1-0.15, 100 times smaller than the
% signals, so it is sure that the ERP response is not the fluctuation, but
% the real response, the next thing is to figure out what does it mean.
%% 4: trial average spectrum
[faxis,Sxxm]=trialAveSpectrum(moveEpoch1,0.001);
figure
subplot(2,2,1)
plot(faxis, Sxxm)
hold on
[faxis,Sxxm]=trialAveSpectrum(unmoveEpoch1,0.001);
plot(faxis, Sxxm)
xlim([0 80]);				%... in select frequency range,
% ylim([0 400])				%... in select power range,
xlabel('Frequency [Hz]')	%... with axes labelled.
ylabel('Power [mV^2/Hz]')
legend('move', 'unmove');
title(['power spectrum TST ',ch1])
legend boxoff
hold off

[faxis,Sxxm]=trialAveSpectrum(moveEpoch2,0.001);
% figure
subplot(2,2,2)
plot(faxis, Sxxm)
hold on
[faxis,Sxxm]=trialAveSpectrum(unmoveEpoch2,0.001);
plot(faxis, Sxxm)
xlim([0 80]);				%... in select frequency range,
% ylim([0 400])				%... in select power range,
xlabel('Frequency [Hz]')	%... with axes labelled.
ylabel('Power [mV^2/Hz]')
legend('move', 'unmove');
title(['power spectrum TST ',ch2])
legend boxoff
hold off

[faxis,Sxxm]=trialAveSpectrum(moveEpoch3,0.001);
% figure
subplot(2,2,3)
plot(faxis, Sxxm)
hold on
[faxis,Sxxm]=trialAveSpectrum(unmoveEpoch3,0.001);
plot(faxis, Sxxm)
xlim([0 80]);				%... in select frequency range,
% ylim([0 400])				%... in select power range,
xlabel('Frequency [Hz]')	%... with axes labelled.
ylabel('Power [mV^2/Hz]')
legend('move', 'unmove');
title(['power spectrum TST ',ch3])
legend boxoff
hold off

[faxis,Sxxm]=trialAveSpectrum(moveEpoch4,0.001);
% figure
subplot(2,2,4)
plot(faxis, Sxxm)
hold on
[faxis,Sxxm]=trialAveSpectrum(unmoveEpoch4,0.001);
plot(faxis, Sxxm)
xlim([0 80]);				%... in select frequency range,
% ylim([0 400])				%... in select power range,
xlabel('Frequency [Hz]')	%... with axes labelled.
ylabel('Power [mV^2/Hz]')
legend('move', 'unmove');
title(['power spectrum TST ',ch4])
legend boxoff
hold off
%% 5: plot Spectrum in dB scale
[faxis,Sxxm]=trialAveSpectrum(moveEpoch1,0.001);
figure
subplot(2,2,1)
plot(faxis, 10*log10(Sxxm/max(Sxxm)))
hold on
[faxis,Sxxm]=trialAveSpectrum(unmoveEpoch1,0.001);
plot(faxis, 10*log10(Sxxm/max(Sxxm)))
xlim([0 80]);				%... in select frequency range,
ylim([-30 0])				%... in select power range,
xlabel('Frequency [Hz]')	%... with axes labelled.
ylabel('Power [dB]')
legend('move', 'unmove');
title(['power spectrum TST ',ch1])
legend boxoff
hold off

[faxis,Sxxm]=trialAveSpectrum(moveEpoch2,0.001);
% figure
subplot(2,2,2)
plot(faxis, 10*log10(Sxxm/max(Sxxm)))
hold on
[faxis,Sxxm]=trialAveSpectrum(unmoveEpoch2,0.001);
plot(faxis, 10*log10(Sxxm/max(Sxxm)))
xlim([0 80]);				%... in select frequency range,
ylim([-30 0])				%... in select power range,
xlabel('Frequency [Hz]')	%... with axes labelled.
ylabel('Power [dB]')
legend('move', 'unmove');
title(['power spectrum TST ',ch2])
legend boxoff
hold off

[faxis,Sxxm]=trialAveSpectrum(moveEpoch3,0.001);
% figure
subplot(2,2,3)
plot(faxis, 10*log10(Sxxm/max(Sxxm)))
hold on
[faxis,Sxxm]=trialAveSpectrum(unmoveEpoch3,0.001);
plot(faxis, 10*log10(Sxxm/max(Sxxm)))
xlim([0 80]);				%... in select frequency range,
ylim([-30 0])				%... in select power range,
xlabel('Frequency [Hz]')	%... with axes labelled.
ylabel('Power [dB]')
legend('move', 'unmove');
title(['power spectrum TST ',ch3])
legend boxoff
hold off

[faxis,Sxxm]=trialAveSpectrum(moveEpoch4,0.001);
% figure
subplot(2,2,4)
plot(faxis, 10*log10(Sxxm/max(Sxxm)))
hold on
[faxis,Sxxm]=trialAveSpectrum(unmoveEpoch4,0.001);
plot(faxis, 10*log10(Sxxm/max(Sxxm)))
xlim([0 80]);				%... in select frequency range,
ylim([-30 0])				%... in select power range,
xlabel('Frequency [Hz]')	%... with axes labelled.
ylabel('Power [dB]')
legend('move', 'unmove');
title(['power spectrum TST ',ch4])
legend boxoff
hold off

%% 6: plot spectrum in dB scale and log-log scale
[faxis,Sxxm]=trialAveSpectrum(moveEpoch1,0.001);
figure
subplot(2,2,1)
semilogx(faxis, 10*log10(Sxxm/max(Sxxm)))	%Log-log scale
hold on
[faxis,Sxxm]=trialAveSpectrum(unmoveEpoch1,0.001);
semilogx(faxis, 10*log10(Sxxm/max(Sxxm)))	%Log-log scale
xlim([0 80]);				%... in select frequency range,
ylim([-30 0])				%... in select power range,
xlabel('Frequency [Hz]')	%... with axes labelled.
ylabel('Power [dB]')
legend('move', 'unmove');
title(['power spectrum TST ',ch1])
legend boxoff
hold off

[faxis,Sxxm]=trialAveSpectrum(moveEpoch2,0.001);
% figure
subplot(2,2,2)
semilogx(faxis, 10*log10(Sxxm/max(Sxxm)))	%Log-log scale
hold on
[faxis,Sxxm]=trialAveSpectrum(unmoveEpoch2,0.001);
semilogx(faxis, 10*log10(Sxxm/max(Sxxm)))	%Log-log scale
xlim([0 80]);				%... in select frequency range,
ylim([-30 0])				%... in select power range,
xlabel('Frequency [Hz]')	%... with axes labelled.
ylabel('Power [dB]')
legend('move', 'unmove');
title(['power spectrum TST ',ch2])
legend boxoff
hold off

[faxis,Sxxm]=trialAveSpectrum(moveEpoch3,0.001);
% figure
subplot(2,2,3)
semilogx(faxis, 10*log10(Sxxm/max(Sxxm)))	%Log-log scale
hold on
[faxis,Sxxm]=trialAveSpectrum(unmoveEpoch3,0.001);
semilogx(faxis, 10*log10(Sxxm/max(Sxxm)))	%Log-log scale
xlim([0 80]);				%... in select frequency range,
ylim([-30 0])				%... in select power range,
xlabel('Frequency [Hz]')	%... with axes labelled.
ylabel('Power [dB]')
legend('move', 'unmove');
title(['power spectrum TST ',ch3])
legend boxoff
hold off

[faxis,Sxxm]=trialAveSpectrum(moveEpoch4,0.001);
% figure
subplot(2,2,4)
semilogx(faxis, 10*log10(Sxxm/max(Sxxm)))	%Log-log scale
hold on
[faxis,Sxxm]=trialAveSpectrum(unmoveEpoch4,0.001);
semilogx(faxis, 10*log10(Sxxm/max(Sxxm)))	%Log-log scale
xlim([0 80]);				%... in select frequency range,
ylim([-30 0])				%... in select power range,
xlabel('Frequency [Hz]')	%... with axes labelled.
ylabel('Power [dB]')
legend('move', 'unmove');
title(['power spectrum TST ',ch4])
legend boxoff
hold off

%% 7: powerspectrogram

x = mean(m2um1);
Fs = 1/dt;					%Define the sampling frequency.
interval = round(Fs);		%Specify the interval size.
overlap = round(Fs*0.95);	%Specify the overlap of intervals.
nfft = round(Fs);			%Specify the FFT length.
[~,F,T,P] = spectrogram(x-mean(x),interval,overlap,nfft,Fs);
figure
subplot(4,2,1)
% imagesc(T,F,10*log10(P))	%... and plot it,
imagesc(T,F,P)
colorbar					%... with a colorbar,
axis xy						%... and origin in lower left, 
hold on
plot([5 5],[0 100],'--','color','w','LineWidth',2)
ylim([0 40])				%... set the frequency range,
xlabel('Time [s]');			%... and label axes.
ylabel('Frequency [Hz]')
title(['TST m2um ',ch1])
caxis([0 20]);
hold off

x = mean(m2um2);
Fs = 1/dt;					%Define the sampling frequency.
interval = round(Fs);		%Specify the interval size.
overlap = round(Fs*0.95);	%Specify the overlap of intervals.
nfft = round(Fs);			%Specify the FFT length.
[~,F,T,P] = spectrogram(x-mean(x),interval,overlap,nfft,Fs);
% figure
subplot(4,2,3)
% imagesc(T,F,10*log10(P))	%... and plot it,
imagesc(T,F,P)
colorbar					%... with a colorbar,
axis xy						%... and origin in lower left, 
hold on
plot([5 5],[0 100],'--','color','w','LineWidth',2)
ylim([0 40])				%... set the frequency range,
xlabel('Time [s]');			%... and label axes.
ylabel('Frequency [Hz]')
title(['TST m2um ',ch2])
caxis([0 20]);
hold off

x = mean(m2um3);
Fs = 1/dt;					%Define the sampling frequency.
interval = round(Fs);		%Specify the interval size.
overlap = round(Fs*0.95);	%Specify the overlap of intervals.
nfft = round(Fs);			%Specify the FFT length.
[~,F,T,P] = spectrogram(x-mean(x),interval,overlap,nfft,Fs);
% figure
subplot(4,2,5)
% imagesc(T,F,10*log10(P))	%... and plot it,
imagesc(T,F,P)
colorbar					%... with a colorbar,
axis xy						%... and origin in lower left, 
hold on
plot([5 5],[0 100],'--','color','w','LineWidth',2)
ylim([0 40])				%... set the frequency range,
xlabel('Time [s]');			%... and label axes.
ylabel('Frequency [Hz]')
title(['TST m2um ',ch3])
caxis([0 100]);
hold off

x = mean(m2um4);
Fs = 1/dt;					%Define the sampling frequency.
interval = round(Fs);		%Specify the interval size.
overlap = round(Fs*0.95);	%Specify the overlap of intervals.
nfft = round(Fs);			%Specify the FFT length.
[~,F,T,P] = spectrogram(x-mean(x),interval,overlap,nfft,Fs);
% figure
subplot(4,2,7)
% imagesc(T,F,10*log10(P))	%... and plot it,
imagesc(T,F,P)
colorbar					%... with a colorbar,
axis xy						%... and origin in lower left, 
hold on
plot([5 5],[0 100],'--','color','w','LineWidth',2)
ylim([0 40])				%... set the frequency range,
xlabel('Time [s]');			%... and label axes.
ylabel('Frequency [Hz]')
title(['TST m2um ',ch4])
caxis([0 40]);
hold off

x = mean(um2m1);
Fs = 1/dt;					%Define the sampling frequency.
interval = round(Fs);		%Specify the interval size.
overlap = round(Fs*0.95);	%Specify the overlap of intervals.
nfft = round(Fs);			%Specify the FFT length.
[~,F,T,P] = spectrogram(x-mean(x),interval,overlap,nfft,Fs);
% figure
subplot(4,2,2)
% imagesc(T,F,10*log10(P))	%... and plot it,
imagesc(T,F,P)
colorbar					%... with a colorbar,
axis xy						%... and origin in lower left, 
hold on
plot([5 5],[0 100],'--','color','w','LineWidth',2)
ylim([0 40])				%... set the frequency range,
xlabel('Time [s]');			%... and label axes.
ylabel('Frequency [Hz]')
title(['TST um2m ',ch1])
caxis([0 20]);
hold off

x = mean(um2m2);
Fs = 1/dt;					%Define the sampling frequency.
interval = round(Fs);		%Specify the interval size.
overlap = round(Fs*0.95);	%Specify the overlap of intervals.
nfft = round(Fs);			%Specify the FFT length.
[~,F,T,P] = spectrogram(x-mean(x),interval,overlap,nfft,Fs);
% figure
subplot(4,2,4)
% imagesc(T,F,10*log10(P))	%... and plot it,
imagesc(T,F,P)
colorbar					%... with a colorbar,
axis xy						%... and origin in lower left, 
hold on
plot([5 5],[0 100],'--','color','w','LineWidth',2)
ylim([0 40])				%... set the frequency range,
xlabel('Time [s]');			%... and label axes.
ylabel('Frequency [Hz]')
title(['TST um2m ',ch2])
caxis([0 20]);
hold off

x = mean(um2m3);
Fs = 1/dt;					%Define the sampling frequency.
interval = round(Fs);		%Specify the interval size.
overlap = round(Fs*0.95);	%Specify the overlap of intervals.
nfft = round(Fs);			%Specify the FFT length.
[~,F,T,P] = spectrogram(x-mean(x),interval,overlap,nfft,Fs);
% figure
subplot(4,2,6)
% imagesc(T,F,10*log10(P))	%... and plot it,
imagesc(T,F,P)
colorbar					%... with a colorbar,
axis xy						%... and origin in lower left, 
hold on
plot([5 5],[0 100],'--','color','w','LineWidth',2)
ylim([0 40])				%... set the frequency range,
xlabel('Time [s]');			%... and label axes.
ylabel('Frequency [Hz]')
title(['TST um2m ',ch3])
caxis([0 100]);
hold off

x = mean(um2m4);
Fs = 1/dt;					%Define the sampling frequency.
interval = round(Fs);		%Specify the interval size.
overlap = round(Fs*0.95);	%Specify the overlap of intervals.
nfft = round(Fs);			%Specify the FFT length.
[S,F,T,P] = spectrogram(x-mean(x),interval,overlap,nfft,Fs);
% figure
subplot(4,2,8)
% imagesc(T,F,10*log10(P))	%... and plot it,
imagesc(T,F,P)
colorbar					%... with a colorbar,
axis xy						%... and origin in lower left, 
hold on
plot([5 5],[0 100],'--','color','w','LineWidth',2)
ylim([0 40])				%... set the frequency range,
xlabel('Time [s]');			%... and label axes.
ylabel('Frequency [Hz]')
title(['TST um2m ',ch4])
caxis([0 40]);
hold off

%% 8: correlation and cohr
% compute the cohr between two electrodes.

[faxis,cohr] = cohrFunc(moveEpoch1,moveEpoch2,1000);
figure
subplot(2,3,1)
plot(faxis, cohr);		%Plot coherence vs frequency,
hold on
[faxis,cohr] = cohrFunc(unmoveEpoch1,unmoveEpoch2,1000);
plot(faxis, cohr)
xlim([0 50])			%... in chosen frequency range,
ylim([0 1])
xlabel('Frequency [Hz]')%... with axes labelled.
ylabel('Coherence')
title(['TST Cohr ',ch1,'-',ch2])
legend('move','unmove')
legend boxoff
hold off

[faxis,cohr] = cohrFunc(moveEpoch1,moveEpoch3,1000);
subplot(2,3,2)
plot(faxis, cohr);		%Plot coherence vs frequency,
hold on
[faxis,cohr] = cohrFunc(unmoveEpoch1,unmoveEpoch3,1000);
plot(faxis, cohr)
xlim([0 50])			%... in chosen frequency range,
ylim([0 1])
xlabel('Frequency [Hz]')%... with axes labelled.
ylabel('Coherence')
title(['TST Cohr ',ch1,'-',ch3])
legend('move','unmove')
legend boxoff
hold off

[faxis,cohr] = cohrFunc(moveEpoch1,moveEpoch4,1000);
subplot(2,3,3)
plot(faxis, cohr);		%Plot coherence vs frequency,
hold on
[faxis,cohr] = cohrFunc(unmoveEpoch1,unmoveEpoch4,1000);
plot(faxis, cohr)
xlim([0 50])			%... in chosen frequency range,
ylim([0 1])
xlabel('Frequency [Hz]')%... with axes labelled.
ylabel('Coherence')
title(['TST Cohr ',ch1,'-',ch4])
legend('move','unmove')
legend boxoff
hold off

[faxis,cohr] = cohrFunc(moveEpoch2,moveEpoch3,1000);
subplot(2,3,4)
plot(faxis, cohr);		%Plot coherence vs frequency,
hold on
[faxis,cohr] = cohrFunc(unmoveEpoch2,unmoveEpoch3,1000);
plot(faxis, cohr)
xlim([0 50])			%... in chosen frequency range,
ylim([0 1])
xlabel('Frequency [Hz]')%... with axes labelled.
ylabel('Coherence')
title(['TST Cohr ',ch2,'-',ch3])
legend('move','unmove')
legend boxoff
hold off

[faxis,cohr] = cohrFunc(moveEpoch2,moveEpoch4,1000);
subplot(2,3,5)
plot(faxis, cohr);		%Plot coherence vs frequency,
hold on
[faxis,cohr] = cohrFunc(unmoveEpoch2,unmoveEpoch4,1000);
plot(faxis, cohr)
xlim([0 50])			%... in chosen frequency range,
ylim([0 1])
xlabel('Frequency [Hz]')%... with axes labelled.
ylabel('Coherence')
title(['TST Cohr ',ch2,'-',ch4])
legend('move','unmove')
legend boxoff
hold off

[faxis,cohr] = cohrFunc(moveEpoch3,moveEpoch4,1000);
subplot(2,3,6)
plot(faxis, cohr);		%Plot coherence vs frequency,
hold on
[faxis,cohr] = cohrFunc(unmoveEpoch3,unmoveEpoch4,1000);
plot(faxis, cohr)
xlim([0 50])			%... in chosen frequency range,
ylim([0 1])
xlabel('Frequency [Hz]')%... with axes labelled.
ylabel('Coherence')
title(['TST Cohr ',ch3,'-',ch4])
legend('move','unmove')
legend boxoff
hold off

%% 9. difference of cohr: move-unmove
[~,cohrm] = cohrFunc(moveEpoch1,moveEpoch2,1000);
[faxis,cohrum] = cohrFunc(unmoveEpoch1,unmoveEpoch2,1000);
figure
subplot(2,3,1)
plot(faxis, cohrm-cohrum)
hold on
plot([0 100],[0 0],'--','color','k')
xlim([0 100])			%... in chosen frequency range,
% ylim([-0.5 0.5])
xlabel('Frequency [Hz]')%... with axes labelled.
ylabel('difference in LFP Coherence')
title(['TST ',ch1,'-',ch2])
hold off

[~,cohrm] = cohrFunc(moveEpoch1,moveEpoch3,1000);
[faxis,cohrum] = cohrFunc(unmoveEpoch1,unmoveEpoch3,1000);
subplot(2,3,2)
plot(faxis, cohrm-cohrum)
hold on
plot([0 100],[0 0],'--','color','k')
xlim([0 100])			%... in chosen frequency range,
% ylim([-0.5 0.5])
xlabel('Frequency [Hz]')%... with axes labelled.
ylabel('difference in LFP Coherence')
title(['TST ',ch1,'-',ch3])
hold off

[~,cohrm] = cohrFunc(moveEpoch1,moveEpoch4,1000);
[faxis,cohrum] = cohrFunc(unmoveEpoch1,unmoveEpoch4,1000);
subplot(2,3,3)
plot(faxis, cohrm-cohrum)
hold on
plot([0 100],[0 0],'--','color','k')
xlim([0 100])			%... in chosen frequency range,
% ylim([-0.5 0.5])
xlabel('Frequency [Hz]')%... with axes labelled.
ylabel('difference in LFP Coherence')
title(['TST ',ch1,'-',ch4])
hold off

[~,cohrm] = cohrFunc(moveEpoch2,moveEpoch3,1000);
[faxis,cohrum] = cohrFunc(unmoveEpoch2,unmoveEpoch3,1000);
subplot(2,3,4)
plot(faxis, cohrm-cohrum)
hold on
plot([0 100],[0 0],'--','color','k')
xlim([0 100])			%... in chosen frequency range,
% ylim([-0.5 0.5])
xlabel('Frequency [Hz]')%... with axes labelled.
ylabel('difference in LFP Coherence')
title(['TST ',ch2,'-',ch3])
hold off

[~,cohrm] = cohrFunc(moveEpoch2,moveEpoch4,1000);
[faxis,cohrum] = cohrFunc(unmoveEpoch2,unmoveEpoch4,1000);
subplot(2,3,5)
plot(faxis, cohrm-cohrum)
hold on
plot([0 100],[0 0],'--','color','k')
xlim([0 100])			%... in chosen frequency range,
% ylim([-0.5 0.5])
xlabel('Frequency [Hz]')%... with axes labelled.
ylabel('difference in LFP Coherence')
title(['TST ',ch2,'-',ch4])
hold off

[~,cohrm] = cohrFunc(moveEpoch3,moveEpoch4,1000);
[faxis,cohrum] = cohrFunc(unmoveEpoch3,unmoveEpoch4,1000);
subplot(2,3,6)
plot(faxis, cohrm-cohrum)
hold on
plot([0 100],[0 0],'--','color','k')
xlim([0 100])			%... in chosen frequency range,
% ylim([-0.5 0.5])
xlabel('Frequency [Hz]')%... with axes labelled.
ylabel('difference in LFP Coherence')
title(['TST ',ch3,'-',ch4])
hold off
% Should add sem shade around. Later do that.
