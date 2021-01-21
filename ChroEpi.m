% This script is for analyse the long time recording of epileptic mice.
% it will be good for detect any spikes.
% called function: touv, segdata, calbs
%% prepare data, step1 read in and touv transform.

[head, dset] = edfread('20210114_GL8665_KA3.5M_DREADD_day7_sal_night_APCR1_dHCR2_S1R3_S1L4.edf');
dset_1 = touv(dset(1,:));   
dset_2 = touv(dset(2,:));   
dset_3 = touv(dset(3,:));   
dset_4 = touv(dset(4,:));   
dset_5 = dset(5,:);
ch1 = 'APCR';                       %...Define channel name
ch2 = 'dHCR';
ch3 = 'S1R';
ch4 = 'S1L';
% define sampling parameters
Fs = 1000; 
dt = 1/Fs; 

%% prepare data step2: scale time length

D = 8*60*60;                        % Define the length of dataset.
dset_1 = dset_1(1:D*Fs);            % slice the data to 8h-long.
dset_2 = dset_2(1:D*Fs);
dset_3 = dset_3(1:D*Fs);
dset_4 = dset_4(1:D*Fs);
dset_5 = dset_5(1:D*Fs);
dset = [dset_1;dset_2;dset_3;dset_4;dset_5]; % rebuild dset after touv and slice.
t = (dt:dt:(size(dset_1,2))*dt);

%% prepare data step3: add filter (bandpass 1-70, fir1 ord=200).

d = designfilt('bandpassiir','FilterOrder',2, ...
    'HalfPowerFrequency1',1,'HalfPowerFrequency2',70, ...
    'SampleRate',1000,'DesignMethod','butter');
                                    % filter design, butterworth method
                                    % ...[1 70], ord = 2.
% fvtool(d)                           % visualize the filter

dset_1 = filtfilt(d,dset_1);        % apply designed filter to ch1
dset_2 = filtfilt(d,dset_2);        % apply designed filter to ch2
dset_3 = filtfilt(d,dset_3);        % apply designed filter to ch3
dset_4 = filtfilt(d,dset_4);        % apply designed filter to ch4


%% calculate baseline of every channel
t0bs = 14130;   % set the t0 of bs calculate, by eye checking on edfBrowser.
durbs = 100;% set the dur of bs calculate, by eye checking on edfBrowser.

x4bs1 = segdata(dset_1,t0bs,1000,durbs);  %Define bs from dset_1.
bs1 = calbs(x4bs1);
x4bs2 = segdata(dset_2,t0bs,1000,durbs);  %Define bs from dset_2.
bs2 = calbs(x4bs2);
x4bs3 = segdata(dset_3,t0bs,1000,durbs);  %Define bs from dset_1.
bs3 = calbs(x4bs3);
x4bs4= segdata(dset_4,t0bs,1000,durbs);   %Define bs from dset_1.
bs4 = calbs(x4bs4);
%% remove extreme value, use values in abs(baseline) to replace.
% x = dset_1;
% [y,i,xmedian,xsigma] = hampel(x,10000,20);
% n = 1:length(x);
% figure;
% plot(n,x)
% hold on
% plot(n,xmedian-3*xsigma,n,xmedian+3*xsigma)
% plot(find(i),x(i),'sk')
% hold off
% % legend('Original signal','Lower limit','Upper limit','Outliers')
% 
% figure;plot(n,y)

% this method also will introduce ineffective peaks.

%% another way to find the replacement 
% replace the founded x by the datapoint in the prepared data.
% x = dset_1;
% [~,i,~,~] = hampel(x,10000,14);
% outlierloc = find(i);
% figure;subplot(2,1,1)
% plot(dset_1)
% ylim([-800 800])
% hold on
% plot(find(i),x(i),'sk')
% x4replace = segdata(x,0,1000,600);
% for n = 1:length(outlierloc)
%     x(outlierloc(n)) = x4replace(round(rand*length(x4replace)));
% end
% subplot(2,1,2)
% plot(x)
% ylim([-800 800])
% dset_1 = x;

% after practice, I think if I do not replace all the outlier peaks with
% other data, the replaced outliers would generate one more peak than the
% original data. So, I think if for spike numbers computation, it is
% unnecessary to do this step.


%% detect epileptic spikes and build the table.
x = dset_2;        %Define the movement channel.
[pks,locs] = findpeaks(-x,t,'MinPeakProminence',2.5*bs2);   
                                %...Define movement threshold, 3*bs.

figure;
plot(t,x)                       %...plot data.
hold on
plot(locs,-pks,'o')             %...plot epileptic peaks
xlabel('time [s]')
ylim([-800 800])
hold off

% detect the spike interval time(> 1s)
difflocs = diff(locs);  %...calculate time interval of the move point.
% figure;
% plot(difflocs)
% ylabel('spike interval[s]')

%% heatmap 
% mat = zeros(8,9);
mat(6,7) = length(locs);

% mat1_5 = zeros(8,0);
% mat1_5(3,1) = length(locs);