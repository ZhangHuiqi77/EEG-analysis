% this script provide common processing for LFP data
% Called hand-made functions: touv, segdata,shuffleFunc, bootstrappingFunc,
% trialAveSpectrum, spectrogramFunc,cohrFunc

%% step 1.plot the data
% read data- edf
[head, dset] = edfread('20201103_GL8211_TST_dHCL1_S1R2_APCR3_APCL4.edf');
dset_1 = touv(dset(3,:)); % APCR %%%%%%%%%%%%%%%%%%%%%%%%%
dset_2 = touv(dset(4,:)); % APCL %%%%%%%%%%%%%%%%%%%%%%%%
dset_3 = touv(dset(1,:)); % dHCR %%%%%%%%%%%%%%%%%%%%%%%%%
dset_4 = touv(dset(2,:)); % S1L %%%%%%%%%%%%%%%%%%%%%%%%%
dset_5 = dset(5,:);
dset = [dset_1;dset_2;dset_3;dset_4;dset_5];
% define sampling parameters
Fs = 1000; %%%%%%%%%%%%%%%%%%%%%%
% gain = 1; %%%%%%%%%%%%%%%%%%%%%%%
dt = 1/Fs; %%%%%%%%%%%%%%%%%%%%%%
t = (dt:dt:(size(dset_1,2))*dt);

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

% xlimit = [100 200]; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ylimit = [-200 200]; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% subplot(5,1,1)
% plot(t,dset_1,'k')
% % xlabel('Time [s]')
% ylabel('Voltage [\mu V]')
% % xlim([xlimit(1) xlimit(2)]);
% ylim([ylimit(1) ylimit(2)]);
% title('raw LFP - APCR')
% 
% subplot(5,1,2)
% plot(t,dset_2,'k')
% % xlabel('Time [s]')
% ylabel('Voltage [\mu V]')
% % xlim([xlimit(1) xlimit(2)]);
% ylim([ylimit(1) ylimit(2)]);
% title('raw LFP - APCL')
% 
% subplot(5,1,3)
% plot(t,dset_3,'k')
% % xlabel('Time [s]')
% ylabel('Voltage [\mu V]')
% % xlim([xlimit(1) xlimit(2)]);
% ylim([ylimit(1) ylimit(2)]);
% title('raw LFP - dHCR')
% 
% subplot(5,1,4)
% plot(t,dset_4,'k')
% % xlabel('Time [s]')
% ylabel('Voltage [\mu V]')
% % xlim([xlimit(1) xlimit(2)]);
% ylim([ylimit(1) ylimit(2)]);
% title('raw EEG - S1L')
% 
% subplot(5,1,5)
% plot(t,dset_5)
% xlabel('Time [s]')
% % xlim([xlimit(1) xlimit(2)]);
% ylabel('acceleration')

%% step2- spectrum 
x = dset_3; %%%%%%%%%%%%%%%%%%%%%%%%% Define interested Ch
t = (dt:dt:(size(x,2))*dt);
N = length(x);
T = t(end);

x = hann(N).*transpose(x);	%Multiply data by Hanning window.
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

%% step3 - split to bands and link axis
x = dset_3; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = t(2)-t(1);				%Define the sampling interval.
Fs = 1/dt;					%Define the sampling frequency.
fNQ = Fs/2;					%Define the Nyquist frequency.
							%For low frequency interval,
Wn = [5,12]/fNQ;				%...set the passband, Theta
n  = 100;					%...and filter order,
b  = fir1(n,Wn);			%...build bandpass filter.
Vlo = filtfilt(b,1,x);	%...and apply filter.
							%For high frequency interval,
Wn = [30,50]/fNQ;			%...set the passband, low Gamma
n  = 100;					%...and filter order,
b  = fir1(n,Wn);			%...build bandpass filter.
Vhi = filtfilt(b,1,x);	%...and apply filter.


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

%% step 4. extract example data
% pick up t0 by hand this time, but this step will be automatically
x = dset_3;
m2um = zeros(10,10*Fs);     % data of 10s for spectrogram, transform
um2m = zeros(10,10*Fs);     % data of 10s for spectrog, transform
t0_m2um = [93 136 147 174 195 232 261 341 294 323]; % prepared by hand
t0_um2m = [105 125 141 181 208 246 271 353 301 328]; % prepared by hand
t0_m2um = sort(t0_m2um);    % Sort the time point
t0_um2m = sort(t0_um2m);    % Sort the time point

for m = 1:10
    m2um(m,:) = segdata(x,t0_m2um(m),Fs,10);
    um2m(m,:) = segdata(x,t0_um2m(m),Fs,10);
end


um_off = um2m(:, 1:0.5*size(um2m,2)); % prepare data of 5s for spectrum
m_on = um2m(:, 0.5*size(um2m,2)+1:end);
m_off = m2um(:, 1:0.5*size(m2um,2));
um_on = m2um(:, 0.5*size(m2um,2)+1:end);
move = [m_on;m_off];
unmove = [um_on;um_off];

%% segment the data automatically
x = dset_5-mean(dset_5);        %Define the movement channel.
Fs = 1000;                      %Define sampling frequency.
dt = 1/Fs;                      %Define time interval for sampling.
t = (dt:dt:(size(x,2))*dt);     %Define time variable. 

x4bs = segdata(x,202,1000,10);  %Define data segment for baseline calculation by eye.
bs = calbs(x4bs);               %Calculate basiline- p2p difference.
[pks_m,locs_m] = findpeaks(x,t,'MinPeakProminence',3*bs);   %Define movement threshold, 3*bs.

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
xlabel('detected moving time interval[s]')
ylabel('[s]')

thr = 1;    %Define a threshold for the unmove duration.
[pks,locs] = findpeaks(difflocs_m,'MinPeakHeight',thr); 
hold on
plot(locs,pks,'o')
t0um = locs_m(locs);    %start time loc of um.
t0m = locs_m(locs+1);   %start time loc of m.
t0m = [0 t0m];          %add 0 to starting point of m.
t0mm = t0m(1:(length(t0m)-1)); %Define t0m as m-um loop. so rm the last one.
t0mu = t0m(2:end);             %Define t0m as um-m loop.

figure;
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
ind = find(mdur==0);           % find the index of mdur is 0, which is too short and won't be selected in the later analysis.
t0mm(ind) = [];             % remain t0mm that dur >=1.
t0mmu = t0um;               % prepare a new t0um;
t0mmu(ind)=[];              % remain coresponding t0um for t0mm
mdur = mdur(mdur>0);    % remain the mdur >=1.

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
moveEpoch2 = data4epoch(dset_2,t0mm,1000,mdur,1);   % xepoch = data4epoch(x,t0,Fs,dur,sec),Ch1
moveEpoch3 = data4epoch(dset_3,t0mm,1000,mdur,1);   % xepoch = data4epoch(x,t0,Fs,dur,sec),Ch1
moveEpoch4 = data4epoch(dset_4,t0mm,1000,mdur,1);   % xepoch = data4epoch(x,t0,Fs,dur,sec),Ch1


%% prepare data: unmove epoch
umdur = t0mu-t0um;              % calculate the duration of every movement. 
umdur = floor(umdur);            % round the mdur to integers.

figure;plot(t,x)    % ...plot and eye checking time for move start.
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
unmoveEpoch2 = data4epoch(dset_2,t0um,1000,umdur,1);   % xepoch = data4epoch(x,t0,Fs,dur,sec),Ch1
unmoveEpoch3 = data4epoch(dset_3,t0um,1000,umdur,1);   % xepoch = data4epoch(x,t0,Fs,dur,sec),Ch1
unmoveEpoch4 = data4epoch(dset_4,t0um,1000,umdur,1);   % xepoch = data4epoch(x,t0,Fs,dur,sec),Ch1

%% prepare data ( >5s um2m)
t0mm = t0m(1:(length(t0m)-1)); %Define t0m as m-um loop. so rm the last one.
t0mu = t0m(2:end);             %Define t0m as um-m loop.
mdur = t0um-t0mm;              % calculate the duration of every movement. 
umdur = t0mu-t0um;              % calculate the duration of every unmovement. 

nloop = length(umdur);     % m-um loop
m2um = zeros(nloop,10*Fs);
um2m = zeros(nloop,10*Fs);
for m = 1:nloop
   if mdur(m)>5 && umdur(m)>5
       m2um(m,:) = segdata(x,t0um(m)-5,Fs,10);
   end
end
for m = 1:nloop-1
    if umdur(m)>5 && mdur(m+1)>5
        um2m(m,:) = segdata(x,t0mu(m)-5,Fs,10);
    end
end
m2um = m2um(any(m2um,2),any(m2um,1));   % remain nonzero lines.
um2m = um2m(any(um2m,2),any(um2m,1));   % remain nonzero lines.
% !well done, then use the auto prepared data to do the ana.

%% statistic 1: ERP and bootstrapping
% plot ERP for the continous 4 stages
x = m2um;
t = (dt:dt:10);
figure;
hold on
for m = 1:size(x,1)
    plot(t,x,'r')
end
plot(t,mean(x),'LineWidth',2,'color','k')
plot([5 5],[-150 200],'LineWidth',1,'color','w');
plot([0 10],[0 0],'LineWidth',1,'color','w')
title('m2um')
hold off

x = um2m;
t = (dt:dt:10);
figure;
hold on
for m = 1:size(x,1)
    plot(t,x,'r')
end
plot(t,mean(x),'LineWidth',2,'color','k')
plot([5 5],[-150 200],'LineWidth',1,'color','w');
plot([0 10],[0 0],'LineWidth',1,'color','w')
title('um2m')
hold off
% the result of oringinal 10 trials shows a little mess and 
% seems no significance in the m2um data, some significant change in um2m
% data. Anyway, bootstrapping give a more convinced result with more repeat
%% plot bootstrapping um2m with shuffled data
dt = 0.001;
t = (-5:dt:5-dt);

xs = shuffleFunc(um2m,1);                   % generate shuffled data
[mn,ciL,ciU] = bootstrappingFunc(xs,3000);  % apply bootstrapping to shuffled data for 3000 times
figure;
plot(t,ciL,'color',[0.8 0.8 0.8])           %... and plot lower CI,
hold on
plot(t,ciU,'color',[0.8 0.8 0.8])           %... and upper CI.
plot(t, mn, 'LineWidth', 3,'color', [0.5 0.5 0.5])    %... and plot it.

% bootstrapping for target data
[mn,ciL,ciU] = bootstrappingFunc(um2m,3000);  % apply bootstrapping for 3000 times
plot(t,ciL,'b')                     %... and plot lower CI,
plot(t,ciU,'b')                     %... and upper CI.
plot(t, mn, 'LineWidth', 3,'color','k')    %... and plot it.

ylabel('Voltage [\mu V]')       %Label the y-axis as voltage.
xlabel('Time [s]')              %Label the x-axis as time.
title('ERP of um2m with bootstrap confidence intervals and shuffled data')
plot([0 0],[-80 80],'LineWidth',1,'color','w');
plot([-5 5],[0 0],'LineWidth',1,'color','w')
hold off

%% plot bootstrappint m2um with shuffled data
dt = 0.001;
t = (-5:dt:5-dt);
% prepare shuffled data
xs = shuffleFunc(m2um,1);                   % generate shuffled data
[mn,ciL,ciU] = bootstrappingFunc(xs,3000);  % apply bootstrapping to shuffled data for 3000 times
figure;
plot(t,ciL,'color',[0.8 0.8 0.8])           %... and plot lower CI,
hold on
plot(t,ciU,'color',[0.8 0.8 0.8])           %... and upper CI.
plot(t, mn, 'LineWidth', 3,'color', [0.5 0.5 0.5])    %... and plot it.

% bootstrapping for target data
[mn,ciL,ciU] = bootstrappingFunc(m2um,3000);  % apply bootstrapping for 3000 times
plot(t,ciL,'b')                     %... and plot lower CI,
plot(t,ciU,'b')                     %... and upper CI.
plot(t, mn, 'LineWidth', 3,'color','k')    %... and plot it.

ylabel('Voltage [\mu V]')       %Label the y-axis as voltage.
xlabel('Time [s]')              %Label the x-axis as time.
title('ERP of m2um with bootstrap confidence intervals and shuffled data')
plot([0 0],[-80 80],'LineWidth',1,'color','w');
plot([-5 5],[0 0],'LineWidth',1,'color','w')
hold off


%% trial average spectrum
[faxis,Sxxm]=trialAveSpectrum(m_on,0.001);
figure
plot(faxis, Sxxm)
hold on
[faxis,Sxxm]=trialAveSpectrum(m_off,0.001);
plot(faxis, Sxxm)

[faxis,Sxxm]=trialAveSpectrum(um_on,0.001);
plot(faxis, Sxxm)

[faxis,Sxxm]=trialAveSpectrum(um_off,0.001);
plot(faxis, Sxxm)
% plot(faxis, 10*log10(Sxx))	%Plot power in decibels vs frequency,
xlim([0 30]);				%... in select frequency range,
ylim([0 400])				%... in select power range,
xlabel('Frequency [Hz]')	%... with axes labelled.
ylabel('Power [ mV^2/Hz]')
legend('move on', 'move off', 'unmove on', 'unmove off');
legend boxoff
hold off
%% plot Sxx in dB scale
[faxis,Sxxm]=trialAveSpectrum(m_on,0.001); % calculate Sxx
figure
plot(faxis, 10*log10(Sxxm/max(Sxxm)))		%Plot power in decibels.
hold on

[faxis,Sxxm]=trialAveSpectrum(m_off,0.001); % calculate Sxx
plot(faxis, 10*log10(Sxxm/max(Sxxm)))		%Plot power in decibels.

[faxis,Sxxm]=trialAveSpectrum(um_on,0.001); % calculate Sxx
plot(faxis, 10*log10(Sxxm/max(Sxxm)))		%Plot power in decibels.

[faxis,Sxxm]=trialAveSpectrum(um_off,0.001); % calculate Sxx
plot(faxis, 10*log10(Sxxm/max(Sxxm)))		%Plot power in decibels.

xlim([0 30])							%Select frequency range.
ylim([-60 0])							%Select decibel range.
xlabel('Frequency [Hz]')				%Label axes.
ylabel('Power [dB]')
legend('move on', 'move off', 'unmove on', 'unmove off');
legend boxoff
hold off

%% plot Sxx in dB scale and log-log scale
figure;
[faxis,Sxxm]=trialAveSpectrum(m_on,0.001); % calculate Sxx
semilogx(faxis, 10*log10(Sxxm/max(Sxxm)))	%Log-log scale
hold on

[faxis,Sxxm]=trialAveSpectrum(m_off,0.001); % calculate Sxx
semilogx(faxis, 10*log10(Sxxm/max(Sxxm)))	%Log-log scale

[faxis,Sxxm]=trialAveSpectrum(um_on,0.001); % calculate Sxx
semilogx(faxis, 10*log10(Sxxm/max(Sxxm)))	%Log-log scale

[faxis,Sxxm]=trialAveSpectrum(um_off,0.001); % calculate Sxx
semilogx(faxis, 10*log10(Sxxm/max(Sxxm)))	%Log-log scale

xlim([df 100])							%Select frequency range.
ylim([-60 0])							%Select decibel range.
xlabel('Frequency [Hz]')				%Label axes.
ylabel('Power [dB]')
legend('move on', 'move off', 'unmove on', 'unmove off');
legend boxoff
hold off

%% powerspectrogram for um2m
spectrogramFunc(um2m,0.001)
hold on
plot([5 5],[0 100],'--','color','w','LineWidth',2)
ylim([0 70])				%... set the frequency range,
xlabel('Time [s]');			%... and label axes.
ylabel('Frequency [Hz]')
title('TST - unmove to move')
colormapeditor
hold off
%% powerspectrogram for um2m
% x = mean(m2um);                %Calculate ERP
spectrogramFunc(um2m(1,:),0.001)
hold on
plot([5 5],[0 100],'--','color','w','LineWidth',2)
ylim([0 70])				%... set the frequency range,
xlabel('Time [s]');			%... and label axes.
ylabel('Frequency [Hz]')
title('TST - move to unmove')
colormapeditor
hold off
% if I image a single line vector in um2m or m2um, i cannot see any
% significant frequecy change.
%% correlation and cohr
% prepare data for 1s as an epoch
moveEpoch1 = data4epoch(dset_1,t0_m2um,1000,4,1);   % xepoch = data4epoch(x,t0,Fs,dur,sec),Ch1
moveEpoch2 = data4epoch(dset_2,t0_m2um,1000,4,1);   % Ch2
moveEpoch3 = data4epoch(dset_3,t0_m2um,1000,4,1);   % Ch3
moveEpoch4 = data4epoch(dset_4,t0_m2um,1000,4,1);   % Ch4

unmoveEpoch1 = data4epoch(dset_1,t0_um2m,1000,4,1);   % prepare data epoch, Ch1
unmoveEpoch2 = data4epoch(dset_2,t0_um2m,1000,4,1);   % prepare data epoch, Ch1
unmoveEpoch3 = data4epoch(dset_3,t0_um2m,1000,4,1);   % prepare data epoch, Ch1
unmoveEpoch4 = data4epoch(dset_4,t0_um2m,1000,4,1);   % prepare data epoch, Ch1

%% compute the cohr between two electrodes.

% E1 = move3;              %Define Electrode 1
% E2 = move4;            %Define Electrode 2
% K = size(E1,1);			%Define the number of trials.
% N = size(E1,2);			%Define the number of indices per trial.
% dt = t(2)-t(1);			%Define the sampling interval.
% t = (dt:dt:1);          %Define the timeline according to epoch.
% T  = t(end); 			%Define the duration of data.
% 
% Sxx = zeros(K,N);		%Create variables to save the spectra,
% Syy = zeros(K,N);
% Sxy = zeros(K,N);
% for k=1:K				%... and compute spectra for each trial.
%     x=E1(k,:)-mean(E1(k,:));
%     y=E2(k,:)-mean(E2(k,:));
%     Sxx(k,:) = 2*dt^2/T * (fft(x) .* conj(fft(x)));
%     Syy(k,:) = 2*dt^2/T * (fft(y) .* conj(fft(y)));
%     Sxy(k,:) = 2*dt^2/T * (fft(x) .* conj(fft(y)));
% end
% 
% Sxx = Sxx(:,1:N/2+1);	%Ignore negative frequencies.
% Syy = Syy(:,1:N/2+1);
% Sxy = Sxy(:,1:N/2+1);
% 
% Sxx = mean(Sxx,1);		%Average the spectra across trials.
% Syy = mean(Syy,1);
% Sxy = mean(Sxy,1);		%... and compute the coherence.
% cohr = abs(Sxy) ./ (sqrt(Sxx) .* sqrt(Syy));
% 
% df = 1/max(T);			%Determine the frequency resolution.
% fNQ = 1/dt/2;			%Determine the Nyquist frequency,
% faxis = (0:df:fNQ);		%... and construct frequency axis.
% 1:APCR; 2:APCL; 3:dHCR; 4:S1L
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
title('TST Cohr APCR-APCL')
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
title('TST Cohr APCR-dHCR')
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
title('TST Cohr APCR-S1L')
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
title('TST Cohr APCL-dHCR')
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
title('TST Cohr APCL-S1L')
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
title('TST Cohr dHCR-S1L')
legend('move','unmove')
legend boxoff
hold off

%% difference of cohr: move-unmove
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
title('TST APCR-APCL')
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
title('TST APCR-dHCR')
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
title('TST APCR-S1L')
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
title('TST APCL-dHCR')
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
title('TST APCL-S1L')
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
title('TST dHCR-S1L')
hold off
% Should add sem shade around. Later do that.


