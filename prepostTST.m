% plot the data
[head, dset] = edfread('20201103_GL8214_preTST_box2_APCR1_LSR2_dHCR3_S1L4.edf');
dset_1 = touv(dset(1,:));   % transfer data to uV.
dset_2 = touv(dset(2,:));   
dset_3 = touv(dset(3,:));   
dset_4 = touv(dset(4,:));   
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

figure('color','w');plot(t,dset_5)
x4bs = segdata(dset_5,193,1000,10);  %Define data for baseline cal by eye.


% plot all 4 channels
figure('color','w');
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
ax1.XLim = [0 30]; 
linkaxes([ax1,ax2,ax3,ax4],'y')
ax1.YLim = [-300 300];  

%% prepare data for pre-TST
% I think I only need to reconstruct it.
x = dset_1-mean(dset_1);
n = t(end)/1;
preTSTepoch1 = reshape(x,[Fs*1,n])';

x = dset_2-mean(dset_2);
n = t(end)/1;
preTSTepoch2 = reshape(x,[Fs*1,n])';

x = dset_3-mean(dset_3);
n = t(end)/1;
preTSTepoch3 = reshape(x,[Fs*1,n])';

x = dset_4-mean(dset_4);
n = t(end)/1;
preTSTepoch4 = reshape(x,[Fs*1,n])';

%% trial-average spectrum
[faxis,Sxxm]=trialAveSpectrum(preTSTepoch1,0.001);
figure('color','w')
subplot(2,2,1)
plot(faxis, Sxxm)
hold on
[faxis,Sxxm]=trialAveSpectrum(postTSTepoch1,0.001);
plot(faxis, Sxxm)
xlim([0 80]);				%... in select frequency range,
% ylim([0 400])				%... in select power range,
xlabel('Frequency [Hz]')	%... with axes labelled.
ylabel('Power [mV^2/Hz]')
legend('pre', 'post');
title(['power spectrum TST ',ch1])
legend boxoff
hold off

[faxis,Sxxm]=trialAveSpectrum(preTSTepoch2,0.001);
% figure
subplot(2,2,2)
plot(faxis, Sxxm)
hold on
[faxis,Sxxm]=trialAveSpectrum(postTSTepoch2,0.001);
plot(faxis, Sxxm)
xlim([0 80]);				%... in select frequency range,
% ylim([0 400])				%... in select power range,
xlabel('Frequency [Hz]')	%... with axes labelled.
ylabel('Power [mV^2/Hz]')
legend('pre', 'post');
title(['power spectrum TST ',ch2])
legend boxoff
hold off

[faxis,Sxxm]=trialAveSpectrum(preTSTepoch3,0.001);
% figure
subplot(2,2,3)
plot(faxis, Sxxm)
hold on
[faxis,Sxxm]=trialAveSpectrum(postTSTepoch3,0.001);
plot(faxis, Sxxm)
xlim([0 80]);				%... in select frequency range,
% ylim([0 400])				%... in select power range,
xlabel('Frequency [Hz]')	%... with axes labelled.
ylabel('Power [mV^2/Hz]')
legend('pre', 'post');
title(['power spectrum TST ',ch3])
legend boxoff
hold off

[faxis,Sxxm]=trialAveSpectrum(preTSTepoch4,0.001);
% figure
subplot(2,2,4)
plot(faxis, Sxxm)
hold on
[faxis,Sxxm]=trialAveSpectrum(postTSTepoch4,0.001);
plot(faxis, Sxxm)
xlim([0 80]);				%... in select frequency range,
% ylim([0 400])				%... in select power range,
xlabel('Frequency [Hz]')	%... with axes labelled.
ylabel('Power [mV^2/Hz]')
legend('pre', 'post');
title(['power spectrum preTST ',ch4])
legend boxoff
hold off
%% coherence
% compute the cohr between two electrodes.

[faxis,cohr] = cohrFunc(preTSTepoch1,preTSTepoch2,1000);
figure('color','w')
subplot(2,3,1)
plot(faxis, cohr);		%Plot coherence vs frequency,
hold on
[faxis,cohr] = cohrFunc(postTSTepoch1,postTSTepoch2,1000);
plot(faxis, cohr)
xlim([0 50])			%... in chosen frequency range,
ylim([0 1])
xlabel('Frequency [Hz]')%... with axes labelled.
ylabel('Coherence')
title(['TST Cohr ',ch1,'-',ch2])
legend('pre','post')
legend boxoff
hold off

[faxis,cohr] = cohrFunc(preTSTepoch1,preTSTepoch3,1000);
subplot(2,3,2)
plot(faxis, cohr);		%Plot coherence vs frequency,
hold on
[faxis,cohr] = cohrFunc(postTSTepoch1,postTSTepoch3,1000);
plot(faxis, cohr)
xlim([0 50])			%... in chosen frequency range,
ylim([0 1])
xlabel('Frequency [Hz]')%... with axes labelled.
ylabel('Coherence')
title(['TST Cohr ',ch1,'-',ch3])
legend('pre','post')
legend boxoff
hold off

[faxis,cohr] = cohrFunc(preTSTepoch1,preTSTepoch4,1000);
subplot(2,3,3)
plot(faxis, cohr);		%Plot coherence vs frequency,
hold on
[faxis,cohr] = cohrFunc(postTSTepoch1,postTSTepoch4,1000);
plot(faxis, cohr)
xlim([0 50])			%... in chosen frequency range,
ylim([0 1])
xlabel('Frequency [Hz]')%... with axes labelled.
ylabel('Coherence')
title(['TST Cohr ',ch1,'-',ch4])
legend('pre','post')
legend boxoff
hold off

[faxis,cohr] = cohrFunc(preTSTepoch2,preTSTepoch3,1000);
subplot(2,3,4)
plot(faxis, cohr);		%Plot coherence vs frequency,
hold on
[faxis,cohr] = cohrFunc(postTSTepoch2,postTSTepoch3,1000);
plot(faxis, cohr)
xlim([0 50])			%... in chosen frequency range,
ylim([0 1])
xlabel('Frequency [Hz]')%... with axes labelled.
ylabel('Coherence')
title(['TST Cohr ',ch2,'-',ch3])
legend('pre','post')
legend boxoff
hold off

[faxis,cohr] = cohrFunc(preTSTepoch2,preTSTepoch4,1000);
subplot(2,3,5)
plot(faxis, cohr);		%Plot coherence vs frequency,
hold on
[faxis,cohr] = cohrFunc(postTSTepoch2,postTSTepoch4,1000);
plot(faxis, cohr)
xlim([0 50])			%... in chosen frequency range,
ylim([0 1])
xlabel('Frequency [Hz]')%... with axes labelled.
ylabel('Coherence')
title(['TST Cohr ',ch2,'-',ch4])
legend('pre','post')
legend boxoff
hold off

[faxis,cohr] = cohrFunc(preTSTepoch3,preTSTepoch4,1000);
subplot(2,3,6)
plot(faxis, cohr);		%Plot coherence vs frequency,
hold on
[faxis,cohr] = cohrFunc(postTSTepoch3,postTSTepoch4,1000);
plot(faxis, cohr)
xlim([0 50])			%... in chosen frequency range,
ylim([0 1])
xlabel('Frequency [Hz]')%... with axes labelled.
ylabel('Coherence')
title(['TST Cohr ',ch3,'-',ch4])
legend('pre','post')
legend boxoff
hold off
