5%% step 1.plot the data
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
gain = 1; %%%%%%%%%%%%%%%%%%%%%%%
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
ax1.XLim = [0 400]; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% prepare shuffled x (x can be any stage, here m2um or um2m)
x = um2m; 
ns = 1;                   %define number of shuffle for each trial as 1, so we can generate a dataframe the same as original data
xs = zeros(ns*size(x,1),size(x,2));			%Vector to hold shuffled x results
for m = 1:size(x,1)
    for n = 1:ns				%For each surrogate,
        xs(n+ns*(m-1),:) = x(m,randperm(size(x,2)));	%Resample by y-index
    end
end

um_off = um2m(:, 1:0.5*size(um2m,2)); % prepare data of 5s for spectrum
m_on = um2m(:, 0.5*size(um2m,2)+1:end);
m_off = m2um(:, 1:0.5*size(m2um,2));
um_on = m2um(:, 0.5*size(m2um,2)+1:end);
move = [m_on;m_off];
unmove = [um_on;um_off];

%% segment the data automatically
x = dset_5-mean(dset_5);
x4bs = segdata(x,200,1000,12);
bs = calbs(x4bs);
[pks,locs] = findpeaks(x,'MinPeakProminence',2*bs);
% figure;
% plot(t,x)
% hold on
% plot(locs/Fs,pks,'o')
difflocs = diff(locs);
% figure;
% plot(difflocs)
% Condition 1) choose the unmove time larger than 2s
% Condition 2) choose the move time larger than 2s
locs_um = zeros(100,2);
locs_m = zeros(100,2);
n = 1;
for m = 1:length(difflocs)
    if difflocs(m) > 2000 && difflocs(m+1) < 2000 % unmove on
        locs_um(n,1) = locs(m); % diff(m) = locs(m+1)-locs(m),so locs(m) is locs_um_on
        locs_um(n,2) = locs(m+1); % as the former, locs(m+1) is the um_off
        locs_m(n,1) = locs_um(n-1,2);
        n = n+1;
    end
end
locs_um_on = nonzeros(locs_um_on)';
for m = 1:(length(locs_um_on)-1)
    if locs_um_on(m+1) - locs_um_on(m) < 1000
        locs_um_on(m+1) = [];
    end
end

    

% 
% locs_m_on = zeros(1,100);
% 
% for m = 1: (length(locs)-1)
%     if locs(m+1) - locs(m) > 2000
%         locs_m_on(m) = locs(m+1);
%     end
% end
% locs_m_on = nonzeros(locs_m_on);
% 
% for m = 1:length(locs_m_on)
%      
%     locs(m+1)-locs(m)
figure;
plot(t,x,locs_m_on/Fs,x(locs_m_on),'o')



%% statistic 1: ERP and bootstrapping
% plot ERP for the continous 4 stages
x = [mean(m_off) mean(um_on)]; % transform from move to unmove
% x = [mean(um_off) mean(m_on)];
t = (dt:dt:2);
figure;plot(t,x,'LineWidth',2)
hold on
for m = 1:size(m_off,1)
    plot(t,[m_off(m,:) um_on(m,:)],'r')
end
% the result shows it is mess and seems no significance in the data
%% plot shuffled x with black and gray
x = xs; % define the data to be bootstrapping
ntrials = size(x,1);
dt = 0.001;
t = (-5:dt:5-dt);
ERP0 = zeros(3000,size(x,2));	%Create empty ERP variable.
for k=1:3000						%For each resampling,
    i=randsample(ntrials,ntrials,1);%... choose the trials,
    EEG0 = x(i,:);				%... create resampled EEG,
    ERP0(k,:) = mean(EEG0,1);		%... save resampled ERP.
end

sERP0=sort(ERP0);               %Sort each column of resampled ERP.
ciL  =sERP0(0.025*size(ERP0,1),:);  %Determine lower CI.
ciU  =sERP0(0.975*size(ERP0,1),:);  %Determine upper CI.

mn = mean(x,1);             %Determine ERP for condition A,
figure;
plot(t,ciL,'color',[0.8 0.8 0.8])      %... and plot lower CI,
hold on
plot(t,ciU,'color',[0.8 0.8 0.8])            %... and upper CI.
plot(t, mn, 'LineWidth', 3,'color', [0.5 0.5 0.5])    %... and plot it.

% from the result, it seems that shuffled signal almost cover the original
% signal, it makes sense that shuffled technique is not good for ERP, but
% it maybe very useful for phase-amplitude analysis.

%% bootstrapping for target data and shuffle
x = um2m; % define the data to be bootstrapping
ntrials = size(x,1);
dt = 0.001;
t = (-5:dt:5-dt);
ERP0 = zeros(3000,size(x,2));	%Create empty ERP variable.
for k=1:3000						%For each resampling,
    i=randsample(ntrials,ntrials,1);%... choose the trials,
    EEG0 = x(i,:);				%... create resampled EEG,
    ERP0(k,:) = mean(EEG0,1);		%... save resampled ERP.
end

sERP0=sort(ERP0);               %Sort each column of resampled ERP.
ciL  =sERP0(0.025*size(ERP0,1),:);  %Determine lower CI.
ciU  =sERP0(0.975*size(ERP0,1),:);  %Determine upper CI.

mn = mean(x,1);             %Determine ERP for condition A,
% hold on                         %Freeze the current plot, 
plot(t,ciL,'b')                     %... and plot lower CI,
plot(t,ciU,'b')                     %... and upper CI.
plot(t, mn, 'LineWidth', 3,'color','k')    %... and plot it.

ylabel('Voltage [\mu V]')       %Label the y-axis as voltage.
xlabel('Time [s]')              %Label the x-axis as time.
title('ERP of um2m with bootstrap confidence intervals and shuffled data')
plot([0 0],[-80 80],'LineWidth',2,'color','w');
plot([-5 5],[0 0],'LineWidth',2,'color','w')
hold off


%% trial average spectrum

E1 = unmove;   
K = size(E1,1);				%Define the number of trials.
N = size(E1,2);				%Define the number of time indices.
dt = 0.001;
t = (dt:dt:5);
T  = t(end);				%Define the duration of data.

Sxx = zeros(K,N);		%Create variable to store each spectrum.
for k=1:K					%For each trial,
    x = E1(k,:);			%... get the data,
    xf  = fft(x-mean(x));	%... compute Fourier transform,
    Sxx(k,:) = 2*dt^2/T *(xf.*conj(xf));%... compute spectrum.
end
Sxx = Sxx(:,1:N/2+1);		%Ignore negative frequencies.
Sxx = mean(Sxx,1);			%Average spectra over trials.

df = 1/max(T);				%Define frequency resolution,
fNQ = 1/dt/2;				%... and Nyquist frequency.
faxis = (0:df:fNQ);			%... to construct frequency axis.

% figure
plot(faxis, Sxx)
% plot(faxis, 10*log10(Sxx))	%Plot power in decibels vs frequency,
xlim([0 100]);				%... in select frequency range,
ylim([0 400])				%... in select power range,
xlabel('Frequency [Hz]')	%... with axes labelled.
ylabel('Power [ mV^2/Hz]')
%% plot in dB scale
x = mean(m_off);
t = (dt:dt:10);
N = length(x);
T = N*dt;
xf = fft(x-mean(x));			%Compute Fourier transform of x.
Sxx = 2*dt^2/T * (xf.*conj(xf));	%Compute power spectrum.
Sxx = Sxx(1:length(x)/2+1);		%Ignore negative frequencies.

df = 1/max(T);					%Determine frequency resolution.
fNQ = 1/ dt / 2;				%Determine Nyquist frequency.
faxis = (0:df:fNQ);				%Construct frequency axis.

% figure
plot(faxis, 10*log10(Sxx/max(Sxx)))		%Plot power in decibels.
xlim([0 100])							%Select frequency range.
ylim([-60 0])							%Select decibel range.
xlabel('Frequency [Hz]')				%Label axes.
ylabel('Power [dB]')

% semilogx(faxis, 10*log10(Sxx/max(Sxx)))	%Log-log scale
% xlim([df 100])							%Select frequency range.
% % ylim([-60 0])							%Select decibel range.
% xlabel('Frequency [Hz]')				%Label axes.
% ylabel('Power [dB]')

%% powerspectrogram
x = m2um; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = mean(x);
Fs = 1/dt;					%Define the sampling frequency.
interval = round(Fs);		%Specify the interval size.
overlap = round(Fs*0.95);	%Specify the overlap of intervals.
nfft = round(Fs);			%Specify the FFT length.

%Compute the spectrogram,
[S,F,T,P] = spectrogram(x-mean(x),interval,overlap,nfft,Fs);
figure
% imagesc(T,F,10*log10(P))	%... and plot it,
imagesc(T,F,P)
colorbar					%... with a colorbar,
axis xy						%... and origin in lower left, 
ylim([0 70])				%... set the frequency range,
xlabel('Time [s]');			%... and label axes.
ylabel('Frequency [Hz]')
title('TST - move to unmove')
colormapeditor

%% correlation and cohr
% prepare the data for cohr, need to include 4Ch in move and unmove state
x = dset_3;
m2um = zeros(10,10*Fs);
um2m = zeros(10,10*Fs);
t0_m2um = [93 136 147 174 195 232 261 341 294 323]; % prepared by hand
t0_um2m = [105 125 141 181 208 246 271 353 301 328]; % prepared by hand
t0_m2um = sort(t0_m2um);
t0_um2m = sort(t0_um2m);

% try 50% move or unmove time point first
t0_move = t0_m2um;  
t0_unmove = t0_um2m;
% hold on
% for k = 1:length(t0_unmove)
%     plot([t0_unmove(k) t0_unmove(k)],[-100 100],'r')
% end

%% prepare data for 1s as an epoch
x = dset_4;     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sec4epoch = 1;      % Define 1s for every epoch in the later corr&cohr ana.
move = zeros(40,sec4epoch*Fs);       % we have 10 variable in t0, and every t0 has include 4s
unmove = zeros(40,sec4epoch*Fs); 
% the former description is because I pick up t0 by hand and every epoch
% include 5s move and 5 sec unmove. Here, for guarantee, we pick 4s from
% each epoch. So totally 40 windows.

for m = 1:length(t0_move)
    toseg = segdata(x,t0_move(m),Fs,4);    % extract move data from x
    move(((m-1)*4+1):(m*4),:) = reshape(toseg,[1000,4])'; % reshape the data to 1s in a window
    toseg = segdata(x,t0_unmove(m),Fs,4);  % extract unmove data from x
    unmove(((m-1)*4+1):(m*4),:) = reshape(toseg,[1000,4])';
end

% % prepare shuffled x (x can be any stage, here m2um or um2m)
% x = move; 
% ns = 1;                   %define number of shuffle for each trial as 1, so we can generate a dataframe the same as original data
% xs = zeros(ns*size(x,1),size(x,2));			%Vector to hold shuffled x results
% for m = 1:size(x,1)
%     for n = 1:ns				%For each surrogate,
%         xs(n+ns*(m-1),:) = x(m,randperm(size(x,2)));	%Resample by y-index
%     end
% end
move4 = move;       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
unmove4 = unmove;       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% compute the cohr between two electrodes.
% 1:APCR; 2:APCL; 3:dHCR; 4:S1L
E1 = move3;              %Define Electrode 1
E2 = move4;            %Define Electrode 2
K = size(E1,1);			%Define the number of trials.
N = size(E1,2);			%Define the number of indices per trial.
dt = t(2)-t(1);			%Define the sampling interval.
t = (dt:dt:1);          %Define the timeline according to epoch.
T  = t(end); 			%Define the duration of data.

Sxx = zeros(K,N);		%Create variables to save the spectra,
Syy = zeros(K,N);
Sxy = zeros(K,N);
for k=1:K				%... and compute spectra for each trial.
    x=E1(k,:)-mean(E1(k,:));
    y=E2(k,:)-mean(E2(k,:));
    Sxx(k,:) = 2*dt^2/T * (fft(x) .* conj(fft(x)));
    Syy(k,:) = 2*dt^2/T * (fft(y) .* conj(fft(y)));
    Sxy(k,:) = 2*dt^2/T * (fft(x) .* conj(fft(y)));
end

Sxx = Sxx(:,1:N/2+1);	%Ignore negative frequencies.
Syy = Syy(:,1:N/2+1);
Sxy = Sxy(:,1:N/2+1);

Sxx = mean(Sxx,1);		%Average the spectra across trials.
Syy = mean(Syy,1);
Sxy = mean(Sxy,1);		%... and compute the coherence.
cohr = abs(Sxy) ./ (sqrt(Sxx) .* sqrt(Syy));

df = 1/max(T);			%Determine the frequency resolution.
fNQ = 1/dt/2;			%Determine the Nyquist frequency,
faxis = (0:df:fNQ);		%... and construct frequency axis.

figure
plot(faxis, cohr);		%Plot coherence vs frequency,
xlim([0 50])			%... in chosen frequency range,
ylim([0 1])
xlabel('Frequency [Hz]')%... with axes labelled.
ylabel('Coherence')
title('Cohr dHCR-S1L TST move')

%% step 