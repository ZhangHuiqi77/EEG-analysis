%% find extreme values
x = dset_1;
[y,i,xmedian,xsigma] = hampel(x,10000,20);
n = 1:length(x);
figure;
plot(n,x)
hold on
plot(n,xmedian-3*xsigma,n,xmedian+3*xsigma)
plot(find(i),x(i),'sk')
hold off
legend('Original signal','Lower limit','Upper limit','Outliers')
figure;plot(n,y)        % this is the replaced data
%% filter design test
% FOR filter test
x = segdata(dset_2,5060,1000,360);
% dt = t(2)-t(1);				%Define the sampling interval.
Fs = 1/dt;					%Define the sampling frequency.

% Wn = [1,70]/fNQ;		    %...set the passband, 1-70
% n  = 100;					%...and filter order,


d1 = designfilt('bandpassiir','FilterOrder',2, ...
    'HalfPowerFrequency1',1,'HalfPowerFrequency2',70, ...
    'SampleRate',1000,'DesignMethod','butter');
% fvtool(d1)
xd1 = filtfilt(d1,x);

% d2 = designfilt('bandpassfir','FilterOrder',2, ...
%          'CutoffFrequency1',1,'CutoffFrequency2',70, ...
%          'SampleRate',1000,'DesignMethod','butter');

d2 = designfilt('bandpassfir', 'FilterOrder', 2, 'CutoffFrequency1', 1, 'CutoffFrequency2', 70, 'SampleRate', 1000, 'DesignMethod', 'window');

% fvtool(d2)
xd2 = filtfilt(d2,x);

% d3 = designfilt('highpassfir','StopbandFrequency',0.25, ...
%          'PassbandFrequency',0.35,'PassbandRipple',0.5, ...
%          'StopbandAttenuation',65,'DesignMethod','kaiserwin');

figure;
subplot(3,1,1)
plot(x)
ylim([-600 200])

subplot(3,1,2)
plot(xd1)
ylim([-600 200])

subplot(3,1,3)
plot(xd2)
ylim([-600 200])

% compare the filter effect of iir and fir filter, I think iir is the most
% effective to remove the low frequency drift and ord = 2 is enough. But it
% will generate some big waves around 5000 data point, I think it is OK
% because all group will have that. 

%% how to creat a function which have optional input variables
% here is an example, if no content input, it will return "hello Matlab",
% if there is some content as an input content, the content will be filled

function [] = speakFunc(name,content)
    if ~exist('content','var')
        content = 'hello MATLAB';
    end
    fprintf(1,'%s speak: %s\n',name,content);
end
%% practice cell and struct

fileID = fopen("D:/Download/Chrome/student",'r','n','UTF-8');
students = {}; % creat a students cell
n = 1; 
studentchar = fgetl(fileID); % read the 1st line
studentcell = strsplit(studentchar,','); % split student info from char to cell
student = struct(); % prepare a empty struct and studentcell features will be in
while isa(studentcell,'cell')
    student.name = studentcell{1};
    student.ID = studentcell{2};
    student.home = studentcell{3};
    student.score = str2double(studentcell{4});
    students(n) = {student};
    n = n+1;
    studentchar = fgetl(fileID);
    if isa(studentchar,'char') % add this iteration, because we know the end of fgetl is -1, which can not be split
        studentcell = strsplit(studentchar,',');
    else
        studentcell = studentchar;
    end
end
%% plot seizure
% define sampling parameters
Fs = 1000; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gain = 1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = 1/Fs; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %read data- edf
% % [head_pre, dset_pre] = edfread('20200524_A8_pilo220_APCR1_APCL2_BLAR3_S1L4.edf');
% % [head_pilo, dset_pilo] = edfread('20200524_A8_pilo220_APCR1_APCL2_BLAR3_S1L4.edf');
[head, dset1] = edfread('20200516_A8_PV_ChETA_bs_rec_APCR1_APCL2_BLAR3_S1L4.edf');
[head, dset2] = edfread('20200524_A8_pilo220_APCR1_APCL2_BLAR3_S1L4.edf');
dset = [dset1 dset2];


dset_1 = dset(1,:); %dHC
dset_2 = dset(2,:); %APCR
dset_3 = dset(3,:); %vHCR
dset_4 = dset(4,:); %BLAR

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


t = [0:dt:(size(dset_1,2)-1)*dt];
N = length(dset_1);
T = N*dt;
df = 1/max(T);
fNQ = 1/dt/2;

figure
hold on
plot(t,dset_1,'k')
plot(t,dset_2-5000,'k')
plot(t,dset_3-5000*2,'k')
plot(t,dset_4-5000*3,'k')
% ylim([-10000,1000])
title('pilo') %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold off

%% calculate power spectrum
x = dset_1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = [0:dt:(size(x,2)-1)*dt];
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
% ylim([0 5000])
xlabel("Frequency [Hz]")
ylabel('Power [ \muV^2/Hz]')
title('spectrum M1L SRS')%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
x = dset_4;
y = double(x);
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
xlim([0 2700])
ylim([0 100]) % ... set the freq range,
xlabel('Time [s]');
ylabel('Frequency [Hz]')
title('APCL pilo') %%%%%%%%%%%%%%%%%


%% clustered bar plot
x = ['pre','light','post'];
pw_alpah_mn = mean(pw_alpha,1);
pw_beta_mn = mean(pw_beta,1);
pw_delta_mn = mean(pw_delta,1);
pw_gamma_mn = mean(pw_gamma,1);
pw_theta_mn = mean(pw_theta,1);
a = pw_delta
b = [pw_delta',pw_theta',pw_alpha',pw_beta',pw_gamma'];
b = double(b);
% define spectrogram parameters
interval = round(Fs); % Specify the interval size.
overlap = round(Fs*0.95); % Specify the overlap of intervals.
nfft = round(Fs); % Specify the FFT length
% have to chanllenge this kind of average

[S,F,T,P] = spectrogram(dset_1,interval, overlap, nfft, Fs);
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
ylim([0 200]) % ... set the freq range,
xlabel('Time [s]');
ylabel('Frequency [Hz]')
title('ChETA Spectrogram- average 6 times under 20mw') %%%%%%%%%%%%%%%%%

hold on
plot([40 40],[0 150],'Color','r')
plot([41 41],[0 150],'Color','r')
plot([40.5 40.5],[0 150],'Color','r')
plot([70.5 70.5],[0 150],'Color','r')
plot(dset_sti)


Fs = 1000;
gain = 10000;
dt = 1/Fs;

hinfo = hdf5info('20200331_A83_KA_3MonPost_M1L1_M1R2_0007');
dset = hdf5read(hinfo.GroupHierarchy.Groups(2).Datasets(1));
dset = dset';
dset_1 = double(dset(1,:)); % M1 L
% dset_1 = [dset_1(600/dt:630/dt) dset_1];
dset_2 = double(dset(2,:)); % M1 R

t = [0:dt:(size(dset_1,2)-1)*dt]; % creat t as x-axis
N = length(dset_1);
T = N*dt;
df = 1/max(T);
fNQ = 1/dt/2;

figure
plot(t,dset_1)
hold on
plot(t,dset_2)
title("A83")
xlabel('Time [s]')
xlim([0 600])
ylim([-50000 50000])

%% read edf file
[head data] = edfread('20200417_A011_WT_piloSE_DHC1_APC2_VHC3_BLA4.edf');
dset_1 = data(1,:); % dHC
dset_2 = data(2,:); % APC
dset_3 = data(3,:); % vHC
dset_4 = data(4,:); % BLA

Fs = 1000;
dt = 1/Fs;

t = [0:dt:(size(dset_1,2)-1)*dt]; % creat t as x-axis
N = length(dset_1);
T = N*dt;
df = 1/max(T);
fNQ = 1/dt/2;

figure
% plot(t,dset_1)
hold on

plot(t,dset_3)
% plot(t,dset_4)
plot(t,dset_2)
title("A011 multi-location recording")
xlabel('Time [s]')
legend('dHC','APC','vHC','BLA')

%% baseline + pilo in APC
[head data] = edfread('20200416_A011_WT_baseline_DHC1_APC2_VHC3_BLA4_1.edf');
baseline = data(2,:); % APC
[head data] = edfread('20200417_A011_WT_piloSE_DHC1_APC2_VHC3_BLA4.edf');
pilo = data(2,:);
dset = [baseline pilo];
Fs = 1000;
dt = 1/Fs;

t = [0:dt:(size(dset,2)-1)*dt]; % creat t as x-axis
t1 = dt*size(baseline,2);
figure
plot(t,dset)
% % plot(t,dset_1)
hold on
% 
% plot(t,dset_3)
% % plot(t,dset_4)
% plot(t,dset_2)
title("20200523_A72_optoSRS_sti_APCR_rec_APCR1_APCL2_BLAR3_S1R4")
xlabel('Time [s]')
plot([t1, t1], [-1000 1000], 'k', 'LineWidth', 2)
%% export files
FolderName = tempdir;   % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = get(FigHandle, 'Name');
  savefig(FigHandle, fullfile(FolderName, 'EEG and analyse post pilo', '.fig'));
end
%% export all plots, use temp filename
 figHandles = findall(0,'Type','figure'); 

 % Loop through figures 2:end
 for i = 1:numel(figHandles)
     fn = tempname("D:\M\output");% creat the output filename
     export_fig(fn, '-tiff', figHandles(i))
 end
 
 
 for k=1:K				%For each trial,
 x = E2(k,:)-mean(E2(k,:));			%...subtract the mean,
 [ac0,lags]=xcorr(x,100,'biased');	%... compute autocovar,
 ac = ac + ac0/K;		%...and add to total, scaled by 1/K.
end
figure
plot(lags*dt,ac)		%Plot autocovar vs lags in time.
xlabel('Lag [s]')		%Label the axes.
ylabel('Autocovariance');


