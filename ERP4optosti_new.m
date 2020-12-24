%% plot Event Related Potential (ERP)
% define sampling parameters
Fs = 1000; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gain = 1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = 1/Fs; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % read data- h5
% hinfo = hdf5info('20200313_A79_CHETA_sti_M1R_rec_M1R1_20HZ_10ms_20mw_30s_180s_times3_0001');
% dset = hdf5read(hinfo.GroupHierarchy.Groups(2).Datasets(1));
% dset = dset';
% dset_1 = double(dset(1,:)); % M1
% dset_1 = [dset_1(600/dt:630/dt) dset_1];
% dset_sti = double(dset(3,:)); % sti tag
% dset_sti = [dset_sti(600/dt:630/dt) dset_sti];
% t = [0:dt:(size(dset_1,2)-1)*dt]; % creat t as x-axis
% N = length(dset_1);
% T = N*dt;
% df = 1/max(T);
% fNQ = 1/dt/2;

%read data- edf
[head, dset] = edfread('20200515_A013_WT_TST_M1R1_M1L2_V1R3_V1L4_1.edf');
dset_1 = dset(1,:); %dHC
dset_2 = dset(2,:); %APCR
dset_3 = dset(3,:); %vHCR
dset_4 = dset(4,:); %BLAR

% hinfo = hdf5info('20200509_A10_PV_ChETA_bs_sti_biAPC_rec_APCL1_APCR2_S1L3_BLAR4_120HZ_5ms_20mw_30s_180s_t6_0001');
% dset_sti = hdf5read(hinfo.GroupHierarchy.Groups(2).Datasets(1));
% dset_sti = double(dset_sti');
% 
% n = (40.5-38.9)*Fs; % correct the lost data in ws at start
% cor = zeros(1,int64(n));
% for i=1:n
%     cor(i)= -300;
% end
% dset_sti = [cor dset_sti];
% 
% n = size(dset_1,2)-size(dset_sti,2); % correct the lost data in ws at last
% cor = zeros(1,n);
% for i=1:n
%     cor(i)= -300;
% end
% dset_sti = [dset_sti cor];

t = (0:dt:(size(dset_1,2)-1)*dt);
N = length(dset_1);
T = N*dt;
df = 1/max(T);
fNQ = 1/dt/2;

% % define opto stimulate parameters
% Fsti = 20; % stimuli frequency [Hz]%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pwidth = 10; % stimuli pulsewidth [ms] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pwidth = pwidth/1000; % convert pwidth to [s]
% dur = 30; % duration [s] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% delay = 180 ;% delay between two stimulations [s]%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cyc = 3; % sti cycles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
hold on
plot(t,dset_1,'k')
plot(t,dset_2-1000,'k')
plot(t,dset_3-1000*2,'k')
plot(t,dset_4-1000*3,'k')
% ylim([-10000,1000])
title('pilo') %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold off

%% locate the Timepoint for each sti and creat sti events
[pks,locs] = findpeaks(dset_sti,'MinPeakHeight',14000,'MinPeakDistance',pwidth*Fs);
ERP0 = zeros(length(locs),pwidth*Fs*4);% prepare a zero matrix for fill in stimuli events
index = int64(locs);
for k = 1:length(locs)
    ERP0(k,:) = dset_1(index(k):index(k)+pwidth*Fs*4-1);
end
ERP1201020 = ERP0/gain*1000; % convert uV to mV, 1:ChETA;0:no-opsin;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate peak-to-peak amptitude
mnERP = mean(ERP1201020,1); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p2pA = max(mnERP(:)) - min(mnERP(:))
p2pA = peak2peak(mnERP);
[M,I] = min(mnERP);
% [M A] = max(mnERP);
latency = I*dt*1000; % [ms] latency = X_lowest-loc_on
a = find(abs(mnERP-(-100))<2);
duration = a*dt*1000;

%% plot ERP for 120hz different power
figure
t2 = (0:dt:pwidth*4-dt);
hold on
plot(t2, mean(ERP01200505,1),'LineWidth',2) % the first 0 means no opsin, then 05-5mw, 10-10mw, 20-20mw
plot(t2, mean(ERP01200510,1),'LineWidth',2)
plot(t2, mean(ERP01200520,1),'LineWidth',2)
plot(t2, mean(ERP1200505,1),'LineWidth',2)
plot(t2, mean(ERP11200510,1),'LineWidth',2)
plot(t2, mean(ERP11200520,1),'LineWidth',2)
xlabel('Time [s]')
ylabel('Voltage [mV]')
ylim([-200 100])
title('LFP trace under 473nm light stimulation in right APC- 120hz 5ms 5/10/20mw') 
legend('no opsin - 5mw','no opsin - 10mw','no opsin - 20mw','ChETA - 5mw','ChETA - 10mw','ChETA - 20mw','Location', 'SouthEast') %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
legend boxoff

%% plot ERP for 20Hz different power
figure
t2 = (0:dt:pwidth*4-dt);
hold on
% for i = 1:k
%     p = plot(t2, ERP0201005(i,:), 'lineWidth',0.5,'color','b');
%     p.Color(4)= 0.05;
% end
% 
% for i = 1:k
%     p = plot(t2, ERP0201010(i,:), 'lineWidth',0.5,'color','g');
%     p.Color(4)= 0.05;
% end
% 
% for i = 1:k
%     p = plot(t2, ERP0201020(i,:), 'lineWidth',0.5,'color','c');
%     p.Color(4)= 0.05;
% end
% 
% for i = 1:k
%     p = plot(t2, ERP1201005(i,:), 'lineWidth',0.5,'color','r');
%     p.Color(4)= 0.05;
% end
% 
% for i = 1:k
%     p = plot(t2, ERP1201010(i,:), 'lineWidth',0.5,'color','k');
%     p.Color(4)= 0.05;
% end
% 
% for i = 1:k
%     p = plot(t2, ERP1201020(i,:), 'lineWidth',0.5,'color','y');
%     p.Color(4)= 0.05;
% end
% plot(t2, mean(ERP0201005,1),'LineWidth',2,'color','b') % the first 0 means no opsin, then 05-5mw, 10-10mw, 20-20mw
% plot(t2, mean(ERP0201010,1),'LineWidth',2,'color','g')
% plot(t2, mean(ERP0201020,1),'LineWidth',2,'color','c')
% plot(t2, mean(ERP1201005,1),'LineWidth',2,'color','r')
% plot(t2, mean(ERP1201010,1),'LineWidth',2,'color','k')
% plot(t2, mean(ERP1201020,1),'LineWidth',2,'color','y')

plot(t2, mean(ERP0201005,1),'LineWidth',2) % the first 0 means no opsin, then 05-5mw, 10-10mw, 20-20mw
plot(t2, mean(ERP0201010,1),'LineWidth',2)
plot(t2, mean(ERP0201020,1),'LineWidth',2)
plot(t2, mean(ERP1201005,1),'LineWidth',2)
plot(t2, mean(ERP1201010,1),'LineWidth',2)
plot(t2, mean(ERP1201020,1),'LineWidth',2)
xlabel('Time [s]')
ylabel('Voltage [mV]')
ylim([-400 100])
title('LFP trace under 473nm light stimulation in right APC- 20hz 10ms 5/10/20mw') 
legend('no opsin - 5mw','no opsin - 10mw','no opsin - 20mw','ChETA - 5mw','ChETA - 10mw','ChETA - 20mw','Location','southeast','NumColumns',1) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot([0 0],[-400 100],'--','Color','k')
% plot([-0.1 -0.1],[-400 100],'--','Color','k')
% plot([-0.05 -0.05],[-400 100],'--','Color','k')
% plot([0.05 0.05],[-400 100],'--','Color','k')
legend boxoff


%% plot ERP for 5HZ
figure
t1 = (-pwidth:dt:pwidth*19-dt);
plot(t1, mean(ERP00501,1),'LineWidth',2,'color','b') % the first 0 means no opsin, then 05-5mw, 10-10mw, 20-20mw
hold on
for i = 1:k
    p = plot(t1, ERP00501(i,:), 'lineWidth',0.5,'color','b');
    p.Color(4)= 0.1;
end
plot(t1, mean(ERP10501,1),'LineWidth',2, 'color','r') % the first 0 means no opsin, then 05-5mw, 10-10mw, 20-20mw
for j = 1:k
    p = plot(t1, ERP10501(j,:), 'lineWidth',0.5,'color','r');
    p.Color(4)= 0.1;
end

plot([0 0],[-2000 1500],'--','Color','k')
xlabel('Time [s]')
ylabel('Voltage [mV]')
xlim([-pwidth pwidth*19])
% ylim([-400 100])
title('LFP trace under 473nm light stimulation in right M1- 5hz 1ms 20mw') 
legend('no opsin','ChETA') %%%%%%%%%%%%%%%%%%%%%
% legend boxoff
hold off

%% plot power spectrum of pre-, light-, post- period
% divide the timeline to pre-, light, post- 30s,calculate power of each period
% see the effect of opto sti to other periods
[pks1,locs1] = findpeaks(dset_sti,t,'MinPeakHeight',14000,'MinPeakDistance',1/Fsti);
loc_on = zeros(1,cyc);
for m=1:cyc
    loc_on(m)= min(locs1);
    locs1(locs1 < loc_on(m)+ dur) = [];
    m = m+1;
end

[pks1,locs1] = findpeaks(dset_sti,t,'MinPeakHeight',14000,'MinPeakDistance',1/Fsti);
loc_off = zeros(1,cyc);
for n=1:cyc
    loc_off(n)= max(locs1);
    locs1(locs1 > loc_off(n)- dur) = [];
    n = n+1;
end
loc_off = sort(loc_off);


% creat every cycle as pre-, light, post- events
x = dset_1; %convert LFP from uV to mV
pre0 = zeros(length(loc_on), dur*Fs);
light0 = zeros(length(loc_on), dur*Fs);
post0 = zeros(length(loc_on), dur*Fs);
for i=1:length(loc_on)
    pre0(i,:) = x(int64((loc_on(i)-dur)/dt+1):int64(loc_on(i)/dt));
    light0(i,:) = x(int64(loc_on(i)/dt+1):int64((loc_on(i)+dur)/dt));
    post0(i,:) = x(int64(loc_off(i)/dt+1):int64((loc_off(i)+dur)/dt));
    i = i+1;
end
% pre = mean(pre0,1);
% light = mean(light0,1);
% post = mean(post0,1);

%% calculate the average power of three events in each period
Sxx_pre = zeros(cyc, dur*Fs/2+1);% R:cycle number; C: is length of Sxx=dur*Fs/2+1
for i = 1:cyc
    x = pre0(i,:);
    t2 = (0:dt:dur-dt);
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

Sxx_light = zeros(cyc, dur*Fs/2+1);
for i = 1:cyc
    x = light0(i,:);
    t2 = (0:dt:dur-dt);
    N = length(x);
    T = N*dt;
    xf = fft(x-mean(x)); % subtract hte mean before compute fft
    xf = double(xf);
    Sxx = 2*dt^2/T*(xf.*conj(xf)); % Compute spectrum
    Sxx = Sxx(1:length(x)/2+1); % Ignore negative freq
    df = 1/max(T); % Determine freq resolution
    fNQ = 1/dt/2; % Determine Nyquist freq
    faxis = (0:df:fNQ); % Construct freq axis
    Sxx_light(i,:) = Sxx;
end

Sxx_post = zeros(cyc, dur*Fs/2+1);
for i = 1:cyc
    x = post0(i,:);
    t2 = (0:dt:dur-dt);
    N = length(x);
    T = N*dt;
    xf = fft(x-mean(x)); % subtract hte mean before compute fft
    xf = double(xf);
    Sxx = 2*dt^2/T*(xf.*conj(xf)); % Compute spectrum
    Sxx = Sxx(1:length(x)/2+1); % Ignore negative freq
    df = 1/max(T); % Determine freq resolution
    fNQ = 1/dt/2; % Determine Nyquist freq
    faxis = (0:df:fNQ); % Construct freq axis
    Sxx_post(i,:) = Sxx;
end

 %% compute bandpower delta, theta, alpha, beta, gamma
% delta 0-4hz
pw_delta = zeros(cyc,3);% Row for repeated times, Col for pre-1 light-2 post-3;
for i = 1:cyc
    X = (df:df:4);
    Y = Sxx_pre(i,df/df:4/df);
    pw_delta(i,1) = trapz(X,Y);
    
    Y = Sxx_light(i,df/df:4/df);
    pw_delta(i,2) = trapz(X,Y);
    
    Y = Sxx_post(i,df/df:4/df);
    pw_delta(i,3) = trapz(X,Y);
end
stage = ["Pre","Light","Post"];
% figure
% boxplot(pw_delta,stage,'BoxStyle','outline','Widths',0.2,'Colors','krb')
% ylabel('Delta band power [mV^2/Hz]')
% % ylim([140000 300000])
% title('no-Opsin sti 20hz 10ms 20mw in M1R-Delta')%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% theta 4-8hz
pw_theta = zeros(cyc,3);% Row for repeated times, Col for pre-1 light-2 post-3;
for i = 1:cyc
    X = (4:df:8);
    Y = Sxx_pre(i,4/df:8/df);
    pw_theta(i,1) = trapz(X,Y);
    
    Y = Sxx_light(i,4/df:8/df);
    pw_theta(i,2) = trapz(X,Y);
    
    Y = Sxx_post(i,4/df:8/df);
    pw_theta(i,3) = trapz(X,Y);
end
stage = ["Pre","Light","Post"];
% figure
% boxplot(pw_theta,stage,'BoxStyle','outline','Widths',0.2,'Colors','krb')
% ylabel('Theta band power [mV^2/Hz]')
% % ylim([140000 300000])
% title('no-Opsin sti 20hz 10ms 20mw in M1R-Theta')%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% alpha 8-13hz
pw_alpha = zeros(cyc,3);% Row for repeated times, Col for pre-1 light-2 post-3;
for i = 1:cyc
    X = (8:df:13);
    Y = Sxx_pre(i,8/df:13/df);
    pw_alpha(i,1) = trapz(X,Y);
    
    Y = Sxx_light(i,8/df:13/df);
    pw_alpha(i,2) = trapz(X,Y);
    
    Y = Sxx_post(i,8/df:13/df);
    pw_alpha(i,3) = trapz(X,Y);
end
stage = ["Pre","Light","Post"];
% figure
% boxplot(pw_alpha,stage,'BoxStyle','outline','Widths',0.2,'Colors','krb')
% ylabel('Alpha band power [mV^2/Hz]')
% % ylim([140000 300000])
% title('no-Opsin sti 20hz 10ms 20mw in M1R-Alpha')%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% beta 13-30hz
pw_beta = zeros(cyc,3);% Row for repeated times, Col for pre-1 light-2 post-3;
for i = 1:cyc
    X = (13:df:30);
    Y = Sxx_pre(i,13/df:30/df);
    pw_beta(i,1) = trapz(X,Y);
    
    Y = Sxx_light(i,13/df:30/df);
    pw_beta(i,2) = trapz(X,Y);
    
    Y = Sxx_post(i,13/df:30/df);
    pw_beta(i,3) = trapz(X,Y);
end
stage = ["Pre","Light","Post"];
% figure
% boxplot(pw_beta,stage,'BoxStyle','outline','Widths',0.2,'Colors','krb')
% ylabel('Beta band power [mV^2/Hz]')
% % ylim([140000 300000])
% title('no-Opsin sti 20hz 10ms 20mw in M1R-Beta')%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Gamma 30-100hz
pw_gamma = zeros(cyc,3);% Row for repeated times, Col for pre-1 light-2 post-3;
for i = 1:cyc
    X = (30:df:100);
    Y = Sxx_pre(i,30/df:100/df);
    pw_gamma(i,1) = trapz(X,Y);
    
    Y = Sxx_light(i,30/df:100/df);
    pw_gamma(i,2) = trapz(X,Y);
    
    Y = Sxx_post(i,30/df:100/df);
    pw_gamma(i,3) = trapz(X,Y);
end
stage = ["Pre","Light","Post"];
% figure
% boxplot(pw_gamma,stage,'BoxStyle','outline','Widths',0.2,'Colors','krb')
% ylabel('Gamma band power [mV^2/Hz]')
% % ylim([140000 300000])
% title('no-Opsin sti 20hz 10ms 20mw in M1R-Gamma')%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% put figures to a subplot
figure

subplot(2,3,1)
boxplot(pw_delta,stage,'BoxStyle','outline','Widths',0.2,'Colors','krb')
ylabel('Delta band power [mV^2/Hz]')
% ylim([140000 300000])
title('ChETA sti 20hz 10ms 20mw in APCR-Delta')%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,3,2)
boxplot(pw_theta,stage,'BoxStyle','outline','Widths',0.2,'Colors','krb')
ylabel('Theta band power [mV^2/Hz]')
% ylim([140000 300000])
title('ChETA sti 20hz 10ms 20mw in APCR-Theta')%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,3,3)
boxplot(pw_alpha,stage,'BoxStyle','outline','Widths',0.2,'Colors','krb')
ylabel('Alpha band power [mV^2/Hz]')
% ylim([140000 300000])
title('ChETA sti 20hz 10ms 20mw in APCR-Alpha')%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,3,4)
boxplot(pw_beta,stage,'BoxStyle','outline','Widths',0.2,'Colors','krb')
ylabel('Beta band power [mV^2/Hz]')
% ylim([140000 300000])
title('ChETA sti 20hz 10ms 20mw in APCR-Beta')%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,3,5)
boxplot(pw_gamma,stage,'BoxStyle','outline','Widths',0.2,'Colors','krb')
ylabel('Gamma band power [mV^2/Hz]')
% ylim([140000 300000])
title('ChETA sti 20hz 10ms 20mw in APCR-Gamma')%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% transfer to graphpad to plot XY
b = [pw_delta',pw_theta',pw_alpha',pw_beta',pw_gamma'];
b = double(b);


%% plot power spectrum
Sxx_pre_mn = mean(Sxx_pre,1);
Sxx_light_mn = mean(Sxx_light,1);
Sxx_post_mn = mean(Sxx_post,1);
figure
% plot(faxis,Sxx_pre_mn, 'k')
% hold on
% plot(faxis, Sxx_light_mn,'r')
% plot(faxis, Sxx_post_mn,'b')
semilogx(faxis, Sxx_pre_mn)
hold on
semilogx(faxis, Sxx_light_mn)
semilogx(faxis, Sxx_post_mn)
xlim([1 100])
% ylim([0 5*10^5])
xlabel("Frequency [Hz]")
ylabel('Power [mV^2/Hz]')
title('PV activated by 120hz optostimuli-A2')%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
legend('pre','light','post')
legend boxoff
hold off

%% plot spectrogram of pre-, light-, post- period
% define spectrogram parameters
interval = round(Fs); % Specify the interval size.
overlap = round(Fs*0.95); % Specify the overlap of intervals.
nfft = round(Fs); % Specify the FFT length
% have to chanllenge this kind of average
pre = mean(pre0,1);
light = mean(light0,1);
post = mean(post0,1);

% Compute the spectrogram
y = [pre light post];
[S,F,T,P] = spectrogram(y-mean(y),interval, overlap, nfft, Fs);
figure
% imagesc(T,F,10*log10(P)) % ... and plot spectrogram
imagesc(T,F,P)
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
ylim([0 130]) % ... set the freq range,
xlabel('Time [s]');
ylabel('Frequency [Hz]')
title('PV activated by 120hz optostimuli-A2') %%%%%%%%%%%%%%%%%