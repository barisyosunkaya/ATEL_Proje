classdef droneRC
%droneRC
%   droneRC is used to create a drone remote controller (RC) object from a
%   .mat file in the MPACT drone RC RF database. 
%
%   droneRC metadata properties:
%       Make
%       Model
%       Index               - Index of the RC in the corresponding folder
%       RawData             - Digitized RF signal
%       ScaleFactor         - Scale factor to convert the digitized signal to volts
%       NumSamples          - Number of samples in the captured signal
%       Duration            - Time span of the captured raw signal
%       Fs                  - Sampling frequency
%       CroppedData         - A portion of the raw signal that includes the signal transient
%                             (to be used instead of raw signal to avoid memory problems)
%       Features            - A structure array of extracted features
%
%   droneRC methods and functions:
%       droneRC             - Create a droneRC object from the database
%       getEnergyTrajectory - Extract normalized energy trajectory from the 
%                             spectrogram of the cropped RF signal
%       findTransientStart  - Find the start point using Higuchi's method
%       plot                - Plot the captured RF signal (or the cropped
%                             portion) and the spectrogram
%
%
% d = droneRC(fileName,Fs) reads the .mat file with the specified name and
% creates a drone RC object d that includes the raw RF signal and the basic
% properties of the captured signal. Fs is the sampling frequency. 
%
% d = droneRC(..., 'OnlyRawData') is the same as d = droneRC(fileName,Fs).
%
% d = droneRC(..., 'OnlyCroppedData') creates a drone RC object d that
% includes only the cropped RF signal and the basic properties.
%
% d = droneRC(..., 'IncludeFeatures') creates a drone RC object d that
% includes the raw and cropped RF signal as well as the extracted features.
% d.Features only includes the start point of the signal transient and can
% be extended to return other features depending on the purpose.
%
% t = getEnergyTrajectory(obj) extracts the normalized energy trajectory
% from the spectrogram of the cropped RF signal. 
%
% s = findTransientStart(obj,W) returns the start point s of the transient from
% the croppped signal using Higuchi's method. W is the size of the window
% for which the fractal dimension is computed.
%
% s = findTransientStart(..., 'Plot') returns the start point s and shows
% it with a line on the cropped RF signal.
%
% plot(obj,'RawData') plots the captured RF signal (Volts) vs time (ms).
% 
% plot(obj,'CroppedData') plots the cropped signal (Volts) vs time (ms).
%
% plot(obj,'Spectrogram') plots the spectrogram of the cropped RF signal.
%
%   Example:
%       % Create a droneRC object
%       fileName = 'MPACT_DroneRC_RF_Database\DJI_Phantom4Pro_1\DJI_Phantom4Pro_1_0001.mat';
%       Fs = 20e9;
%       d = droneRC(fileName,Fs,'IncludeFeatures')
%       W = 50;
%       t = findTransientStart(d,W,'Plot')
% 
% Copyright MPACT Lab, Electrical and Computer Engineering, NC State
% University

    properties
        Make
        Model
        Index
        RawData
        ScaleFactor
        NumSamples % samples
        Duration % s
        Fs = 20e9;% Hz
        CroppedData = []; %samples
        Features = [];
    end
    
    %=====================================================================
    methods
        
        %Constructor
        function obj = droneRC(fileName,Fs,varargin)
            narginchk(0,3)
            
            if nargin == 0
                return
            end
            
            splittedFileName = strsplit(fileName, '\');
            makeModelID = splittedFileName{end};
            splittedMakeModelID = strsplit(makeModelID, '_');
            obj.Make = splittedMakeModelID{1};
            if numel(splittedMakeModelID) == 4
                obj.Model = [splittedMakeModelID{2} '-' splittedMakeModelID{3}];
                obj.Index = splittedMakeModelID{4}(1:end-4);
            else
                obj.Model = splittedMakeModelID{2};
                obj.Index = splittedMakeModelID{3}(1:end-4);
            end
            
            load(fileName,'data')
            
            obj.RawData = double(data.Data-mean(data.Data));
            obj.ScaleFactor = data.YInc;
            crop_start = 2*10^6;
            crop_end = 3*10^6;
            obj.CroppedData = droneRC.cropData(obj.RawData,crop_start,crop_end);
            obj.Fs = Fs;
            obj.NumSamples = data.NumPoints;
            obj.Duration = obj.NumSamples/obj.Fs;
            
            if isempty(varargin) || ( ~isempty(varargin) && strcmp(varargin{1},'OnlyRawData') )
                obj.CroppedData = [];
            elseif ( ~isempty(varargin) && strcmp(varargin{1},'OnlyCroppedData') )
                obj.RawData = [];
            elseif ( ~isempty(varargin) && strcmp(varargin{1},'IncludeFeatures') )
                W = 50;
                obj.Features.TransientStart = findTransientStart(obj,W);
            else
                error('Invalid option!')
            end
            
        end % end constructor 
        
        % Extract normalized energy trajectory from sprectrogram
        function trajectory = getEnergyTrajectory(obj)
            [s,~,~] = spectrogram(double(obj.CroppedData),128,120,256,obj.Fs);
            trajectory = (max(abs(s),[],1))/norm(max(abs(s),[],1));
        end
        
        % Find start point of the RF signal transient
        function transientStart = findTransientStart(obj,W,varargin)
            narginchk(2,3)
            
            fd_arr = [];
            data = double(obj.CroppedData);
            for j = 1:floor(length(data)/W)
                fd_arr(j) = droneRC.Higuchi(data(W*(j-1)+1:W*j));
            end
            
            transientStart = findchangepts(fd_arr)*W;
            
            if ~isempty(varargin) && strcmp(varargin{1},'Plot')
                figure
                plot(obj,'CroppedData')
                hold on
                xline(transientStart/obj.Fs*1e3,'r','LineWidth',1.5)
                xlabel('Time (ms)')
                ylabel('Amplitude (V)')
                title([obj.Make '-' obj.Model '-' obj.Index])
            end
        end
        
        % Plot raw data/cropped data and spectrogram of the cropped data
        function plot(obj,varargin)
            narginchk(2,2)
            
            plotOpts = {'RawData', 'CroppedData', 'Spectrogram'};
            switch varargin{1}
                case plotOpts
                    % Plot raw data
                    if strcmp(varargin{1},plotOpts{1})
                        plot( (1:numel(obj.RawData)) / obj.Fs*1e3, obj.RawData*obj.ScaleFactor )
                        xlabel('Time (ms)')
                        ylabel('Amplitude (V)')
                        title([obj.Make '-' obj.Model '-' obj.Index])
                    % Plot cropped data
                    elseif strcmp(varargin{1},plotOpts{2})
                        plot( (1:numel(obj.CroppedData)) / obj.Fs*1e3, obj.CroppedData*obj.ScaleFactor )
                        xlabel('Time (ms)')
                        ylabel('Amplitude (V)')
                        title([obj.Make '-' obj.Model '-' obj.Index])
                    % plot spectrogram
                    elseif strcmp(varargin{1},plotOpts{3})
                        spectrogram(double(obj.CroppedData),128,120,256,obj.Fs);
                        title([obj.Make '-' obj.Model '-' obj.Index])
                    end
                otherwise
                    error('Invalid feature!');
            end
        end
        
    end       
    
    %=====================================================================
    methods (Static=true,Access=private)
        
        %Crop raw data
        function res = cropData(data,start_point,end_point)
            res = data(start_point:end_point);
        end
        
        % Calculate fractal dimension (FD) using Higuchi's method
        function fd = Higuchi(x)       
            N = length(x);
            kmax = 5; % maximal degree of reduction
            Lmk = zeros(kmax,kmax);
            
            % Get the mean length of curve
            for k = 1:kmax
                for m = 1:k
                    Lmki = 0;
                    for i=1:fix((N-m)/k)
                        Lmki = Lmki+abs(x(m+i*k)-x(m+(i-1)*k));
                    end
                    Ng =        (N-1)/(fix((N-m)/k)*k);
                    Lmk(m,k) =  (Lmki*Ng)/k;
                end
            end
            
            Lk = zeros(1,kmax);
            for k = 1:kmax
                Lk(1,k) = sum(Lmk(1:k,k))/k;
            end
            
            % Calculate the logarithmic values for slope calculation
            lnLk = log(Lk);
            lnk = log(1./[1:kmax]);
            
            % Calculate the slope and assign it to output
            b = polyfit(lnk,lnLk,1);
            fd = b(1);
        end
        
    end
    
end