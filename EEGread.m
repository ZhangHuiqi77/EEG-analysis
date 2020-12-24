function [] = EEGread()
% This function is for plot raw EEG data, acceleration and spectrogram for the selected channels
%   For the Medusa system, we always record four channels and the location
%   can be find in the file name. This function is for us to quickly open
%   and plot the channels we are interested in and wants to go on the later
%   analyse. 

%   If you want to go on with spectrum, spectrogram, power analyse, you
%   need to return at least data from one channel. If you want to go on
%   with coherence analyse, you need ot return two channels.
[fn,fp,index] = uigetfile('*.edf');
    if index == 0
        disp('No file selected')
    else
        [head, dset] = edfread(strcat(fp,fn));
    end
    
    % define data sampling parameter
    
    Fs = 1000; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gain = 1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dt = 1/Fs; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %ask for the channel you want to plot and define the acceleration ch
    prompt = 'which channel you want to plot?';
    ch = input(prompt);
    x = dset(ch,:);
    x = (x*4.5)/(2^23-1)*10^6; % convert the AD value to uV
    acce = dset(5,:);

    % plot data.
    t = [dt:dt:(size(x,2))*dt];
    figure
    subplot(3,1,1)
    plot(t,x,'k')
    xlabel('time [s]')
    ylabel('Voltage [uV]')
    formatSpec = 'Raw signal for channel %d';
    title(sprintf(formatSpec,ch))
    
    % plot Acceleration
    subplot(3,1,2)
    plot(t,acce)
    xlabel('time [s]')
    ylabel('Acceleration')   
    
    % Compute and plot the spectrogram
    interval = round(Fs); % Specify the interval size.
    overlap = round(Fs*0.95); % Specify the overlap of intervals.
    nfft = round(Fs); % Specify the FFT length
    
    y = double(x);
    [S,F,T,P] = spectrogram(y-mean(y),interval, overlap, nfft, Fs);
    subplot(3,1,3)
%     imagesc(T,F,10*log10(P)) % ... and plot spectrogram
    imagesc(T,F,P)
    colormap jet
    colorbar
    set(gca,'colorscale','log')

    axis xy % ... and origin in lower left,
    ylim([0 50])
    xlabel('Time [s]');
    ylabel('Frequency [Hz]')
    formatSpec = 'Spectrogram for channel %d';
    title(sprintf(formatSpec,ch))

end

