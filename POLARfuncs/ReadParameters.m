function [Parameters] = ReadParameters(flnm,FLDR)

% LAST EDITED: 13 May 2021
%==========================================================================
% DESCRIPTION:
%   Reads text file called "PARAMETERS.txt." Information in PARAMETERS.txt 
%   is useful for displaying measurement parameters on Mueller matrix graph,
%   and it tells MATLAB whether to use transmission or reflection mode
%   calculations for the Mueller matrix.
%   Transmission mode is used for the vortex retarders.
%
% REFERENCE:
%   ...
%
% CREATED BY:
%   Matt Goforth, 16 Apr 2021
%==========================================================================

%% OPEN/READ TEXT FILE
%==========================================================================
flpth = fullfile(FLDR,flnm);
F1 = fopen(flpth,'r'); 

CNT = 1; % initialize
while ~feof(F1)
    LINE{CNT} = fgetl(F1);
    CNT = CNT + 1;
end
fclose(F1);

%% POLARIMETER PARAMETERS
%==========================================================================
DATE = strtrim(strsplit(LINE{2},':'));
    Parameters.DATE = DATE{2};
SAMPLE = strtrim(strsplit(LINE{3},':'));
    Parameters.SAMPLE = SAMPLE{2};
MODE = strtrim(strsplit(LINE{4},':'));
    Parameters.MODE = MODE{2};
STEP = strtrim(strsplit(LINE{5},':'));
    Parameters.STEP = STEP{2};
RESOLUTION = strtrim(strsplit(LINE{6},':'));
    Parameters.RESOLUTION = RESOLUTION{2};
LAMBDAo = strtrim(strsplit(LINE{7},':'));
    Parameters.LAMBDAo = LAMBDAo{2};
LAMBDAf = strtrim(strsplit(LINE{8},':'));
    Parameters.LAMBDAf = LAMBDAf{2};
INCIDENT = strtrim(strsplit(LINE{9},':'));
    Parameters.INCIDENT = INCIDENT{2};   
SENSITIVITY = strtrim(strsplit(LINE{10},':'));
    Parameters.SENSITIVITY = SENSITIVITY{2};   
REFIDX = strtrim(strsplit(LINE{11},':'));
    Parameters.REFIDX = REFIDX{2};
end

