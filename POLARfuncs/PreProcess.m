function [Angles,Wavelengths,Intensities,NumWL,NumANG] = PreProcess(InptPATH)

% Last Edited 25 Mar 2021
%==========================================================================
% DESCRIPTION:
%   Separates raw polarimeter data into arrays of angles, wavelengths, and irradiances
% 
% REFERENCE:
%   ...
%
% CREATED BY:
%   Matt Goforth, 09 Sept 2020
%==========================================================================

%% READ IN DATA
%==========================================================================
RawData = dlmread(InptPATH,'\t',0,1); % read raw polarimeter data into array

% Terminate if... 
if min(min(RawData(:,2:end))) < 0
    error('Negative Intensity!'); % ...a negative intensity is found
elseif max(max(RawData(:,2:end))) > 9.9
    error('Lock-In Amp Saturated!');  % ...an intensity above 9.9 is found
end

%% SEPARATE INTO...
%==========================================================================
Wavelengths = round(RawData(:,1),1); % wavelengths measured

                                     % round wavelengths to nearest tenths place;
                                     % for some reason they are slightly off from
                                     % values specified in LabVIEW UI
                                   
Intensities = RawData(:,2:end); % array of measured intensities
[NumWL,NumANG] = size(Intensities); % find number of wavelengths and angles
Angles = 0:6:(6*(NumANG-1)); % angles assumed by R1 during measurement;...
                                ...R2 takes on 5*Angles
end
