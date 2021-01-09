%% step 1.plot the data
% read data- edf for 3rdHS
[head, dset] = edfread('20201227_GL8672_KA3M_DREADD_day1_sal_APCR1_dHCR2_S1R3_S1L4_raw.edf');
dset_1 = dset(1,:); % APCR %%%%%%%%%%%%%%%%%%%%%%%%%
dset_2 = dset(2,:); % dHCR %%%%%%%%%%%%%%%%%%%%%%%%
dset_3 = dset(3,:); % S1R %%%%%%%%%%%%%%%%%%%%%%%%%
dset_4 = dset(4,:); % S1L %%%%%%%%%%%%%%%%%%%%%%%%%
dset_5 = dset(9,:); % volocity, in raw data structure
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
% title('raw LFP - dHCR')

ax3 = nexttile;
plot(t,dset_3,'k')
% title('raw EEG - S1R')

ax4 = nexttile;
plot(t,dset_4,'k')
% title('raw EEG - S1L')
% xlabel('Time [s]')
ylabel('Voltage [\mu V]')

ax5 = nexttile;
plot(t,dset_5)
% title('acceleration')
xlabel('Time [s]')

linkaxes([ax1,ax2,ax3,ax4,ax5],'x')
ax1.XLim = [0 500]; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
linkaxes([ax1,ax2,ax3,ax4],'y')
ax1.YLim = [-300 300];  %%%%%%%%%%%%%%%%%

%% step2- spectrum 
x = dset_2; %%%%%%%%%%%%%%%%%%%%%%%%% Define interested Ch
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
x = dset_2; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = t(2)-t(1);				%Define the sampling interval.
t = (dt:dt:(size(dset_1,2))*dt);
Fs = 1/dt;					%Define the sampling frequency.
fNQ = Fs/2;					%Define the Nyquist frequency.
							%For low frequency interval,
Wn = [5,12]/fNQ;				%...set the passband, Theta
n  = 100;					%...and filter order,
b  = fir1(n,Wn);			%...build bandpass filter.
Vlo = filtfilt(b,1,x);	%...and apply filter.
							%For high frequency interval,
Wn = [30,48]/fNQ;			%...set the passband, low Gamma
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
%% trial average spectrum
Fs = 1000; %%%%%%%%%%%%%%%%%%%%%%
gain = 1; %%%%%%%%%%%%%%%%%%%%%%%
dt = 1/Fs; %%%%%%%%%%%%%%%%%%%%%%
t = (dt:dt:(size(dset_1,2))*dt);
E1 = reshape(dset_2,[1*Fs,t(end)])';   
E1 = E1(1:500,:);
K = size(E1,1);				%Define the number of trials.
N = size(E1,2);				%Define the number of time indices.
dt = 0.001;
t = (dt:dt:1);              %%%%%%%% Define the epoch window
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

figure
plot(faxis, Sxx)
% plot(faxis, 10*log10(Sxx))	%Plot power in decibels vs frequency,
xlim([0 100]);				%... in select frequency range,
ylim([0 400])				%... in select power range,
xlabel('Frequency [Hz]')	%... with axes labelled.
ylabel('Power [ mV^2/Hz]')
%% plot in dB scale
figure
semilogx(faxis, 10*log10(Sxx/max(Sxx)))	%Log-log scale
xlim([df 100])							%Select frequency range.
% ylim([-60 0])							%Select decibel range.
xlabel('Frequency [Hz]')				%Label axes.
ylabel('Power [dB]')

%% powerspectrogram
x = E1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
title('HS-3rd')
colormapeditor
%% plot spectrogram in dB scale
figure
imagesc(T,F,10*log10(P))	%... and plot it,
colorbar					%... with a colorbar,
axis xy						%... and origin in lower left, 
ylim([0 70])				%... set the frequency range,
xlabel('Time [s]');			%... and label axes.
ylabel('Frequency [Hz]')
title('HS-3rd')
colormapeditor
