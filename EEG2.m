%% plot seizure
% define sampling parameters
Fs = 1000; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gain = 1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = 1/Fs; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% read data- edf
[head, dset] = edfread('20200515_A013_WT_TST_M1R1_M1L2_V1R3_V1L4_1.edf');
% last 4min

% [head, dset1] = edfread('20200516_A8_PV_ChETA_bs_rec_APCR1_APCL2_BLAR3_S1L4.edf');
% [head, dset2] = edfread('20200524_A8_pilo220_APCR1_APCL2_BLAR3_S1L4.edf');
% dset = [dset1(:,1:600*Fs) dset2];
dset_1 = dset(1,180*Fs:end); %APCR
dset_2 = dset(2,180*Fs:end); %APCL
dset_3 = dset(3,180*Fs:end); %BLAR
dset_4 = dset(4,180*Fs:end); %S1L

% 
% % read data- h5
% hinfo = hdf5info('20200331_A83_KA_3MonPost_M1L1_M1R2_0003_morning');
% dset = hdf5read(hinfo.GroupHierarchy.Groups(2).Datasets(1));
% dset = dset';
% dset_1 = double(dset(1,:)); % M1
% dset_2 = double(dset(2,:)); % M1
% t = [0:dt:(size(dset_1,2)-1)*dt]; % creat t as x-axis
% N = length(dset_1);
% T = N*dt;
% df = 1/max(T);
% fNQ = 1/dt/2;


t = (0:dt:(size(dset_1,2)-1)*dt);
N = length(dset_1);
T = N*dt;
df = 1/max(T);
fNQ = 1/dt/2;

figure
hold on
plot(t,dset_1,'k')
plot(t,dset_2-1000,'k')
plot(t,dset_3-1000*2,'k')
plot(t,dset_4-1000*3,'k')
% ylim([-10000,1000])
title('pilo') %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold off

%% calculate power spectrum
x = dset_3; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = (0:dt:(size(x,2)-1)*dt);
N = length(x);
T = N*dt;
xf = fft(x-mean(x)); % subtract hte mean before compute fft
xf = double(xf);
Sxx = 2*dt^2/T*(xf.*conj(xf)); % Compute spectrum
Sxx = Sxx(1:length(x)/2+1); % Ignore negative freq
df = 1/max(T); % Determine freq resolution
fNQ = 1/dt/2; % Determine Nyquist freq
faxis = (0:df:fNQ); % Construct freq axis

figure
semilogx(faxis, Sxx, 'k')
xlim([1 100])
ylim([0 5000])
xlabel("Frequency [Hz]")
ylabel('Power [ \muV^2/Hz]')
title('spectrum V1R TST')%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plot spectrum as dB scale
figure
plot(faxis, 10*log10(Sxx/max(Sxx))) % Plot spectrum in decibels-scaled by dividing the max
xlim([0 200]) % Select frequency range
ylim([-100 0]) % Select decibel range
xlabel('Frequency [Hz]') % Label axes
ylabel('Power [dB]')
title('spectrum')

%% plot the spectrogram
Fs = 1/dt; % Define the sampling frequency
interval = round(Fs); % Specify the interval size.
overlap = round(Fs*0.95); % Specify the overlap of intervals.
nfft = round(Fs); % Specify the FFT length

% Compute the spectrogram
x = dset_1;
y = double(x);
[S,F,T,P] = spectrogram(y-mean(y),interval, overlap, nfft, Fs);
figure
imagesc(T,F,10*log10(P)) % ... and plot spectrogram
% imagesc(T,F,P)
colormap jet
colorbar
set(gca,'colorscale','log')
            % cb = colorbar();
            % cb.Ruler.Scale = 'log';
            % cb.Ruler.MinorTick = 'on';
% colormapeditor
% ax = gca;
% mymap = colormap(ax);
% save('MyColormap','mymap') % 10 100000
axis xy % ... and origin in lower left,
% xlim([0 2700])
ylim([0 100]) % ... set the freq range,
xlabel('Time [s]');
ylabel('Frequency [Hz]')
title('M1R TST') %%%%%%%%%%%%%%%%%

%% draw spectrum in one figure
%% calculate mean power spectrum
% during SE, rand to choose ten time point and calculate mean power
%ctrl
r = 300 + 270*rand(10,1);
rt = int64(r*Fs);
dset_ctrl = zeros(10,30*Fs);
for i = 1:10
    dset_ctrl(i,:)= dset_2(rt(i):rt(i)+30*Fs-1);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

Sxx_ctrl = zeros(10, 30*Fs/2+1);% R:repeat number; C: is length of Sxx=dur*Fs/2+1
for i = 1:10
    x = dset_ctrl(i,:);
    t2 = [0:dt:30-dt];
    N = length(x);
    T = N*dt;
    xf = fft(x-mean(x)); % subtract hte mean before compute fft
    xf = double(xf);
    Sxx = 2*dt^2/T*(xf.*conj(xf)); % Compute spectrum
    Sxx = Sxx(1:length(x)/2+1); % Ignore negative freq
    df = 1/max(T); % Determine freq resolution
    fNQ = 1/dt/2; % Determine Nyquist freq
    faxis = (0:df:fNQ); % Construct freq axis
    Sxx_ctrl(i,:) = Sxx;
end
Sxx_ctrl_mn = mean(Sxx_ctrl,1);

% pre-ictal
r = 800 + 170*rand(10,1);
rt = int64(r*Fs);
dset_pre = zeros(10,30*Fs);
for i = 1:10
    dset_pre(i,:)= dset_2(rt(i):rt(i)+30*Fs-1);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

Sxx_pre = zeros(10, 30*Fs/2+1);% R:repeat number; C: is length of Sxx=dur*Fs/2+1
for i = 1:10
    x = dset_pre(i,:);
    t2 = (0:dt:30-dt);
    N = length(x);
    T = N*dt;
    xf = fft(x-mean(x)); % subtract hte mean before compute fft
    xf = double(xf);
    Sxx = 2*dt^2/T*(xf.*conj(xf)); % Compute spectrum
    Sxx = Sxx(1:length(x)/2+1); % Ignore negative freq
    df = 1/max(T); % Determine freq resolution
    fNQ = 1/dt/2; % Determine Nyquist freq
    faxis = (0:df:fNQ); % Construct freq axis
    Sxx_pre(i,:) = Sxx;
end
Sxx_pre_mn = mean(Sxx_pre,1);

% ictal
r = 2250 + 750*rand(10,1);
rt = int64(r*Fs);
dset_pilo = zeros(10,30*Fs);
for i = 1:10
    dset_pilo(i,:)= dset_2(rt(i):rt(i)+30*Fs-1);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

Sxx_pilo = zeros(10, 30*Fs/2+1);% R:repeat number; C: is length of Sxx=dur*Fs/2+1
for i = 1:10
    x = dset_pilo(i,:);
    t2 = (0:dt:30-dt);
    N = length(x);
    T = N*dt;
    xf = fft(x-mean(x)); % subtract hte mean before compute fft
    xf = double(xf);
    Sxx = 2*dt^2/T*(xf.*conj(xf)); % Compute spectrum
    Sxx = Sxx(1:length(x)/2+1); % Ignore negative freq
    df = 1/max(T); % Determine freq resolution
    fNQ = 1/dt/2; % Determine Nyquist freq
    faxis = (0:df:fNQ); % Construct freq axis
    Sxx_pilo(i,:) = Sxx;
end
Sxx_pilo_mn = mean(Sxx_pilo,1);


figure
semilogx(faxis, Sxx_ctrl_mn,'k','LineWidth',2)
hold on
semilogx(faxis, Sxx_pre_mn,'b','LineWidth',2)
semilogx(faxis, Sxx_pilo_mn,'r','LineWidth',2)
% xlim([1 100])
% ylim([0 300000])
xlabel("Frequency [Hz]")
ylabel('Power [ \muV^2/Hz]')
% title('spectrum APCL mean')%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
legend('control','preictal','ictal')
legend boxoff
ax = gca;
ax.FontSize = 16;
% xlim([25 80])
% ylim([0 1600])

% % plot loglog power spectrum
% figure
% loglog(faxis, Sxx_ctrl_mn,'k','LineWidth',2)
% hold on
% loglog(faxis, Sxx_pilo_mn,'b','LineWidth',2)
% % xlim([1 100])
% % ylim([0 300000])
% xlabel("Frequency [Hz]")
% ylabel('Power [ \muV^2/Hz]')
% % title('spectrum APCL mean')%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% legend('preictal','ictal')
% legend boxoff
% ax = gca;
% ax.FontSize = 16;
% xlim([0.1 200])
% % ylim([0 1600])
%%
x = dset_4(1:30*Fs);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t2 = (0:dt:30-dt);
T = 30;
xf = fft(x-mean(x)); % subtract hte mean before compute fft
xf = double(xf);
Sxx = 2*dt^2/T*(xf.*conj(xf)); % Compute spectrum
Sxx = Sxx(1:length(x)/2+1); % Ignore negative freq
df = 1/max(T); % Determine freq resolution
fNQ = 1/dt/2; % Determine Nyquist freq
faxis = (0:df:fNQ); % Construct freq axis
Sxx_pre = Sxx;

x = dset_4(1000*Fs:1030*Fs-1);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t2 = [0:dt:30-dt];
T = 30;
xf = fft(x-mean(x)); % subtract hte mean before compute fft
xf = double(xf);
Sxx = 2*dt^2/T*(xf.*conj(xf)); % Compute spectrum
Sxx = Sxx(1:length(x)/2+1); % Ignore negative freq
df = 1/max(T); % Determine freq resolution
fNQ = 1/dt/2; % Determine Nyquist freq
faxis = (0:df:fNQ); % Construct freq axis
Sxx_pilo = Sxx;

figure
semilogx(faxis, Sxx_pre, 'k','LineWidth',2)
hold on
semilogx(faxis, Sxx_pilo, 'b','LineWidth',2)
xlim([1 100])
% ylim([0 5000])
xlabel("Frequency [Hz]")
ylabel('Power [ \muV^2/Hz]')
title('spectrum S1L')%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
legend('preictal','ictal')
legend boxoff
% set(gca, 'box', 'off')

%% Speed track
% speed = textscan(fileread('TST_speed.csv'),'%f');
fid = fopen('TST_speed.csv');
C = textscan(fid,'%f%f','HeaderLines',8,'Delimiter',',','EndOfLine','\r\n','ReturnOnError',false);
fclose(fid);
[Time,Speed] = C{:}
% Time = [241/7112:241/7112:241];
TIme = [((N-1)/dt)/length(Speed):((N-1)/dt)/length(Speed):(N-1)/dt];% aline the time of speed ana and power ana
figure
plot(Time,Speed,'LineWidth',2,'color','k')