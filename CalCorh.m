% calculate coherence between two electrodes
% when data include multi trials, it is better to use hanning taper
% when data has only 1 trial, it is better to use multi taper
%% plot the two electrodes
[head, dset] = edfread('20200903_GL6834_DREADDexp_pilo290_APCR1_APCL2_S1R3_S1L4.edf');
dset_1 = dset(1,:); 
dset_2 = dset(2,:); 
dset_3 = dset(3,:); 
dset_4 = dset(4,:); 

% define sampling parameters
Fs = 1000; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gain = 1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = 1/Fs;

% define the prepared dataset
x = dset_1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = touv(x); % transfer the value to uV
y = dset_3; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = touv(y);
t = (dt:dt:(size(x,2))*dt);

figure
subplot(2,1,1);
plot(t,x,'k')
title('APCR') %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,1,2);
plot(t,y,'k')
title('S1R') %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% build E1 and E2 in 3 stage, pregs,gs,se
%% manully annotate the stage
t0 = [50 802 832 2000]; %%%%%%%%%%%%%%%%%% annotate the onset time of each stage, baseline, pre_GS, GS, SE
E1 = segdata(x,t0,Fs,30);
E2 = segdata(y,t0,Fs,30);


%%
dt_event = 10;
E1 = zeros(100,dt_event*Fs);
E2 = zeros(100,dt_event*Fs);

t_on = 956;
t_off = t_on+dt_event;
for m = 1:size(E1,1)
    E1(m,:)= x((t_on/dt):(t_off/dt-1));
    E2(m,:)= y((t_on/dt):(t_off/dt-1));
    t_on = t_off+dt;
    t_off = t_on + dt_event;
end



%% calculate cohr - rectangular taper
K = size(E1,1);			%Define the number of trials.
N = size(E1,2);			%Define the number of indices per trial.
t = (dt:dt:dt_event);  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
title("20200903 GL6834 exp")

%% corh by hanning taper instead of rectangular taper
K = size(E1,1);			%Define the number of trials.
N = size(E1,2);			%Define the number of indices per trial.
t = (dt:dt:dt_event);
T  = t(end); 			%Define the duration of data.

Sxx = zeros(K,N);		%Create variables to save the spectra,
Syy = zeros(K,N);
Sxy = zeros(K,N);
for k=1:K				%... and compute spectra for each trial.
    x=E1(k,:)-mean(E1(k,:));
    x=hann(N)'.*x;      %... apply Hanning taper to each trial
    y=E2(k,:)-mean(E2(k,:));
    y=hann(N)'.*y;
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
xlim([0 100])			%... in chosen frequency range,
ylim([0 1])
xlabel('Frequency [Hz]')%... with axes labelled.
ylabel('Coherence hanning taper')
title("20200903 GL6834 ctrl1")

%% cohr multi-taper
corh = zeros(size(E1,1),size(E1,2)/2+1);
for n=1:size(E1,1)
    [C,phi,S12,S1,S2,f]=cohrmt(E1(n,:),E2(n,:),Fs,30);
    corh(n,:) = C';
end

plot(f,C)				%Plot the coherence vs frequency,
ylim([0 1])				%... set the vertical axis,
xlabel('Frequency [Hz]')%... and label the axes.
ylabel('Coherence')
legend("pre GS","GS","SE","baseline")
legend boxoff
title("20200903 GL6834 multitaper for 1 trial")