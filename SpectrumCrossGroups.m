% this is for DREADD experiment comparision
%%
Fs = 1000; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gain = 1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = 1/Fs; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Sxx_ctrl1 = zeros(size(gs_ctrl1,1), 10*Fs/2+1);% R:repeat number; C: is length of Sxx=dur*Fs/2+1
for i = 1:size(gs_ctrl1,1)
    x = gs_ctrl1(i,:);
    t1 = (dt:dt:10);
    N = length(x);
    T = N*dt;
    xf = fft(x-mean(x)); % subtract hte mean before compute fft
    xf = double(xf);
    Sxx = 2*dt^2/T*(xf.*conj(xf)); % Compute spectrum
    Sxx = Sxx(1:length(x)/2+1); % Ignore negative freq
    df = 1/max(T); % Determine freq resolution
    fNQ = 1/dt/2; % Determine Nyquist freq
    faxis = (0:df:fNQ); % Construct freq axis
    Sxx_ctrl1(i,:) = Sxx;
end
Sxx_ctrl1_mn = mean(Sxx_ctrl1,1);


Sxx_ctrl2 = zeros(size(gs_ctrl2,1), 10*Fs/2+1);% R:repeat number; C: is length of Sxx=dur*Fs/2+1
for i = 1:size(gs_ctrl2,1)
    x = gs_ctrl2(i,:);
    t1 = (dt:dt:10);
    N = length(x);
    T = N*dt;
    xf = fft(x-mean(x)); % subtract hte mean before compute fft
    xf = double(xf);
    Sxx = 2*dt^2/T*(xf.*conj(xf)); % Compute spectrum
    Sxx = Sxx(1:length(x)/2+1); % Ignore negative freq
    df = 1/max(T); % Determine freq resolution
    fNQ = 1/dt/2; % Determine Nyquist freq
    faxis = (0:df:fNQ); % Construct freq axis
    Sxx_ctrl2(i,:) = Sxx;
end
Sxx_ctrl2_mn = mean(Sxx_ctrl2,1);


Sxx_exp = zeros(size(gs_exp,1), 10*Fs/2+1);% R:repeat number; C: is length of Sxx=dur*Fs/2+1
for i = 1:size(gs_exp,1)
    x = gs_exp(i,:);
    t1 = (dt:dt:10);
    N = length(x);
    T = N*dt;
    xf = fft(x-mean(x)); % subtract hte mean before compute fft
    xf = double(xf);
    Sxx = 2*dt^2/T*(xf.*conj(xf)); % Compute spectrum
    Sxx = Sxx(1:length(x)/2+1); % Ignore negative freq
    df = 1/max(T); % Determine freq resolution
    fNQ = 1/dt/2; % Determine Nyquist freq
    faxis = (0:df:fNQ); % Construct freq axis
    Sxx_exp(i,:) = Sxx;
end
Sxx_exp_mn = mean(Sxx_exp,1);

figure
% semilogx(faxis, Sxx_ctrl1_mn,'k','LineWidth',2)
% hold on
% semilogx(faxis, Sxx_ctrl2_mn,'b','LineWidth',2)
% semilogx(faxis, Sxx_exp_mn,'r','LineWidth',2)

plot(faxis, Sxx_ctrl1_mn,'k','LineWidth',2)
hold on
plot(faxis, Sxx_ctrl2_mn,'b','LineWidth',2)
plot(faxis, Sxx_exp_mn,'r','LineWidth',2)

xlim([1 100])
% ylim([0 300000])
xlabel("Frequency [Hz]")
ylabel('Power [ \muV^2/Hz]')
% title('spectrum APCL mean')%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
legend('no virus+CNO','hM3Dq+saline','hM3Dq+CNO')
legend boxoff
ax = gca;
ax.FontSize = 16;
% xlim([25 80])
% ylim([0 1600])

%% combine ctrl1 and ctrl2
gs_ctrl = [gs_ctrl1;gs_ctrl2];

Sxx_ctrl = zeros(size(gs_ctrl,1), 10*Fs/2+1);% R:repeat number; C: is length of Sxx=dur*Fs/2+1
for i = 1:size(gs_ctrl,1)
    x = gs_ctrl(i,:);
    t1 = (dt:dt:10);
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

figure
% semilogx(faxis, Sxx_ctrl1_mn,'k','LineWidth',2)
% hold on
% semilogx(faxis, Sxx_ctrl2_mn,'b','LineWidth',2)
% semilogx(faxis, Sxx_exp_mn,'r','LineWidth',2)

plot(faxis, Sxx_ctrl_mn,'k','LineWidth',2)
hold on
plot(faxis, Sxx_exp_mn,'r','LineWidth',2)

xlim([1 100])
% ylim([0 300000])
xlabel("Frequency [Hz]")
ylabel('Power [ \muV^2/Hz]')
% title('spectrum APCL mean')%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
legend('ctrl','hM3Dq+CNO')
legend boxoff
ax = gca;
ax.FontSize = 16;
% xlim([25 80])
% ylim([0 1600])

%% compute bandpower for mean Sxx
pw_delta = zeros(1,3);
X = (df:df:4); % delta 0-4hz
Y = Sxx_ctrl1_mn((df/df):(4/df));
pw_delta(1) = trapz(X,Y); %ctrl1
Y = Sxx_ctrl2_mn((df/df):(4/df));
pw_delta(2) = trapz(X,Y); %ctrl2
Y = Sxx_exp_mn((df/df):(4/df));
pw_delta(3) = trapz(X,Y); %exp

pw_theta = zeros(1,3);
X = (4:df:8);% theta 4-8hz
Y = Sxx_ctrl1_mn((4/df):(8/df));
pw_theta(1) = trapz(X,Y); %ctrl1
Y = Sxx_ctrl2_mn((4/df):(8/df));
pw_theta(2) = trapz(X,Y); %ctrl2
Y = Sxx_exp_mn((4/df):(8/df));
pw_theta(3) = trapz(X,Y); %exp
 
pw_alpha = zeros(1,3);
X = (8:df:13);% alpha 8-13hz
Y = Sxx_ctrl1_mn((8/df):(13/df));
pw_alpha(1) = trapz(X,Y); %ctrl1
Y = Sxx_ctrl2_mn((8/df):(13/df));
pw_alpha(2) = trapz(X,Y); %ctrl2
Y = Sxx_exp_mn((8/df):(13/df));
pw_alpha(3) = trapz(X,Y); %exp

pw_beta = zeros(1,3);
X = (13:df:30);% beta 13-30hz
Y = Sxx_ctrl1_mn((13/df):(30/df));
pw_beta(1) = trapz(X,Y); %ctrl1
Y = Sxx_ctrl2_mn((13/df):(30/df));
pw_beta(2) = trapz(X,Y); %ctrl2
Y = Sxx_exp_mn((13/df):(30/df));
pw_beta(3) = trapz(X,Y); %exp

pw_gamma = zeros(1,3);
X = (30:df:100);% gamma 30-100hz
Y = Sxx_ctrl1_mn((30/df):(100/df));
pw_gamma(1) = trapz(X,Y); %ctrl1
Y = Sxx_ctrl2_mn((30/df):(100/df));
pw_gamma(2) = trapz(X,Y); %ctrl1
Y = Sxx_exp_mn((30/df):(100/df));
pw_gamma(3) = trapz(X,Y); %ctrl1

%%
% delta 0-4hz
pw_delta = zeros(size(Sxx_ctrl1,1),3);
% Row for repeated times(The biggest in the three groups);
% Col for ctrl1 ctrl2 exp;
for m = 1:size(Sxx_ctrl1,1)
    X = (df:df:4);
    Y = Sxx_ctrl1(m,(df/df):(4/df));
    pw_delta(m,1) = trapz(X,Y);
end
for m = 1:size(Sxx_ctrl2,1)
    Y = Sxx_ctrl2(m,(df/df):(4/df));
    pw_delta(m,2) = trapz(X,Y);
end
for m = 1:size(Sxx_exp,1)
    Y = Sxx_exp(m,(df/df):(4/df));
    pw_delta(m,3) = trapz(X,Y);
end

% theta 4-8hz
pw_theta = zeros(size(Sxx_ctrl1,1),3);
% Row for repeated times(The biggest in the three groups);
% Col for ctrl1 ctrl2 exp;
for m = 1:size(Sxx_ctrl1,1)
    X = (4:df:8);
    Y = Sxx_ctrl1(m,(4/df):(8/df));
    pw_theta(m,1) = trapz(X,Y);
end
for m = 1:size(Sxx_ctrl2,1)
    Y = Sxx_ctrl2(m,(4/df):(8/df));
    pw_theta(m,2) = trapz(X,Y);
end
for m = 1:size(Sxx_exp,1)
    Y = Sxx_exp(m,(4/df):(8/df));
    pw_theta(m,3) = trapz(X,Y);
end

% alpha 8-13hz
pw_alpha = zeros(size(Sxx_ctrl1,1),3);
% Row for repeated times(The biggest in the three groups);
% Col for ctrl1 ctrl2 exp;
for m = 1:size(Sxx_ctrl1,1)
    X = (8:df:13);
    Y = Sxx_ctrl1(m,(8/df):(13/df));
    pw_alpha(m,1) = trapz(X,Y);
end
for m = 1:size(Sxx_ctrl2,1)
    Y = Sxx_ctrl2(m,(8/df):(13/df));
    pw_alpha(m,2) = trapz(X,Y);
end
for m = 1:size(Sxx_exp,1)
    Y = Sxx_exp(m,(8/df):(13/df));
    pw_alpha(m,3) = trapz(X,Y);
end

% beta 13-30hz
pw_beta = zeros(size(Sxx_ctrl1,1),3);
% Row for repeated times(The biggest in the three groups);
% Col for ctrl1 ctrl2 exp;
for m = 1:size(Sxx_ctrl1,1)
    X = (13:df:30);
    Y = Sxx_ctrl1(m,(13/df):(30/df));
    pw_beta(m,1) = trapz(X,Y);
end
for m = 1:size(Sxx_ctrl2,1)
    Y = Sxx_ctrl2(m,(13/df):(30/df));
    pw_beta(m,2) = trapz(X,Y);
end
for m = 1:size(Sxx_exp,1)
    Y = Sxx_exp(m,(13/df):(30/df));
    pw_beta(m,3) = trapz(X,Y);
end

%Gamma 30-100hz
pw_gamma = zeros(size(Sxx_ctrl1,1),3);
% Row for repeated times(The biggest in the three groups);
% Col for ctrl1 ctrl2 exp;
for m = 1:size(Sxx_ctrl1,1)
    X = (30:df:100);
    Y = Sxx_ctrl1(m,(30/df):(100/df));
    pw_gamma(m,1) = trapz(X,Y);
end
for m = 1:size(Sxx_ctrl2,1)
    Y = Sxx_ctrl2(m,(30/df):(100/df));
    pw_gamma(m,2) = trapz(X,Y);
end
for m = 1:size(Sxx_exp,1)
    Y = Sxx_exp(m,(30/df):(100/df));
    pw_gamma(m,3) = trapz(X,Y);
end

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