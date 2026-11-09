function droneRC_RFDatabase = createDatabase(databasePath,Fs,NumSignals,opts,varargin)
%createDatabase creates a database of drone remote controllers (RCs)
%
% db = createDatabase(databasePath,Fs,NumSignals,'OnlyRawData') creates a
% database db of drone RC objects that includes only the raw data and the
% basic properties of each drone RC in the subfolders of the specified
% path. Fs is the sampling frequency. NumSignals is the desired number of
% signals from each drone RC folder. 
%
% db = createDatabase(...,'OnlyCroppedData') creates a database db that
% includes only the cropped data and the basic properties to save space.
%
% db = createDatabase(...,'IncludeFeatures') creates a database db that
% contains all the available data. 
%
% db = createDatabase(...,'Table') creates a database db in table format.
%
%   Example:
%       % Create a drone RC database in table format
%       databasePath = 'MPACT_DroneRC_RF_Database\';
%       Fs = 20e9;
%       NumSignals = 10;
%       db = createDatabase(databasePath, Fs, NumSignals, 'IncludeFeatures', 'Table')
%
% Copyright MPACT Lab, Electrical and Computer Engineering, NC State
% University

narginchk(4,5)

folderList = dir(databasePath);
folderList = folderList(3:end);

NumDroneRCs = length(folderList);

droneTable = table;
fields_ = properties(droneRC);

for droneType = 1:NumDroneRCs
    addpath([databasePath folderList(droneType).name '\']);
    fileList = dir([databasePath folderList(droneType).name '\', '*.mat']);
    
    for selectSignal = 1:NumSignals
        clc;
        % Display progress
        fprintf('Creating database...\nController:%d/%d... %%%d\n',droneType,...
            NumDroneRCs,round((droneType-1)/NumDroneRCs*100+ selectSignal/NumSignals/NumDroneRCs*100));
        newEntry = droneRC(fileList(selectSignal).name,Fs,opts);
        
        if ~isempty(varargin) && strcmp(varargin{1}, 'Table')
            
            row_ = cell(1,numel(fields_));
            for p = 1:numel(fields_)
                row_{p} = newEntry.(fields_{p});
            end
            
            droneTable = [droneTable;cell2table(row_)];
            
        elseif isempty(varargin)
            droneRC_RFDatabase{droneType,selectSignal} = newEntry;
            
        else
            error('Invalid input argument!')
        end
        
    end
end

if ~isempty(varargin) && strcmp(varargin{1}, 'Table')
    droneTable.Properties.VariableNames = fields_;
    droneRC_RFDatabase = droneTable;
end

end