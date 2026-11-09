% Create a drone RC RF database in table format with 10 signals from each
% drone make&model
databasePath = uigetdir(); % select the folder "MPACT_DroneRC_RF_Database"
Fs = 20e9; % Hz
NumSignals = 1;

% RFdatabase = createDatabase([databasePath
% '\'],Fs,NumSignals,'IncludeFeatures');
RFdatabase = createDatabase([databasePath '\'],Fs,NumSignals,...
    'IncludeFeatures', 'Table');
head(RFdatabase,7) % display first seven entries

%%
% Create a droneRC object
basepath="C:\Users\ruzga\Documents\GitHub Repositories\ATEL_Proje\matlab files"
%fileName = 'E:\MPACT_DroneRC_RF_Dataset\MPACT_DroneRC_RF_Dataset\DJI_Phantom4Pro_1\DJI_Phantom4Pro_1_0001.mat';
fileName = basepath+"\DJI_Phantom4Pro_1_0001.mat";
% "C:\Users\ruzga\Documents\GitHub Repositories\ATEL_Proje\matlab files\DJI_Phantom4Pro_1_0001.mat"
%denemedenemedenemedeneme
Fs = 20e9;
d = droneRC(fileName,Fs,'IncludeFeatures')
W = 5e3;
t = findTransientStart(d,W,'Plot')

figure(2);
plot(d,'RawData') %plots the captured RF signal (Volts) vs time (ms).
figure(4);
plot(d,'CroppedData') %plots the cropped signal (Volts) vs time (ms).
figure(5);
plot(d,'Spectrogram') %plots the spectrogram of the cropped RF signal.