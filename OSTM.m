% This script is for calculate the sniffing time in OSTM
D = readmatrix('D:\Neu\DATA\videos\Arthur\OSTM_GL8661_1stDLC_resnet_50_OSTMFeb18shuffle1_100000.csv');
Fs = 30;
objx = mean(D(:,11));
objy = mean(D(:,12));

dist = sqrt((D(:,11)-D(:,2)).^2+(D(:,12)-D(:,3)).^2); % distence from nose to the obj.
D(:,14) = dist;
D(:,15) = [diff(dist);0];
dispd = D(2229:2239,:);
disp(dispd)

%step1: set the sniffing threshold
% from my eye check, 70 may be a good threshold
thresh = 50;        % defined by eye
%step2: sniffing conditions: 1)distence<threshold; 2)diff>0
% Condition 1)
FrameNum = length(dist(dist<thresh));
Frame2Sec = FrameNum/Fs;
disp(Frame2Sec);