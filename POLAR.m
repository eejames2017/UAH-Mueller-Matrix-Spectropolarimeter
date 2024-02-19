% Last Edited 07 June 2021
%==========================================================================
% DESCRIPTION: 
%   For 'TOG = 1', specify calibration file and construct system's
%   corrected data reduction matrix 'DRM'. Calibration errors found using 
%   method of least-squares fit. 
%
%   For 'TOG = 0', apply 'DRM' to specified sample data file. Calculate sample's:
%       (1) Raw/uncorrected Mueller matrix (RMM).
%       (2) Calibration-corrected Mueller matrix (MM).
%       (3) Nearest physically realizable Mueller matrix (PMM), found using normalized 'MM'.
%       (4) Deconstructed 'PMM', giving constituent polarization Mueller
%           matrices that together makeup 'PMM'.
%       (5) Polarization properties.
%
% CREATED BY: 
%   Prashant Raman, UAH Physics PhD Student, 2012
% 
% EDITED BY: 
%   Matthew Goforth, UAH Physics Master's Student, 15 Aug 2020
%   
% EDITS: 
%   (09Sep2020) -Added comments.
%               -Grouped all functions into 'POLARfuncs' folder.
%               -Added path to access functions.
%               -Deleted unnecessary files/folders from 'POLAR' folder.
%               -Added for-loop over Mueller matrix plots, where needed. 
%               -Cleaned/organized all functions.
%               -Created 'PreProcess' function to pull irradiances, wavelengths, and angles from raw polarimeter data. 
%               -Altered functions to take in polarimeter irradiances & wavelength values; 
%                   before alteration, each function individually pulled these values from a folder.
%   (22Oct2020) -'lsqcurvefit' function used to fit polarimeter Mueller matrix to
%                   measurements now specifically uses 'levenberg-marquardt' method; 
%                   Sahar's dissertation states this is the fitting method that should be used.
%   (01Dec2020) -Added error warning if an element from the normalized physical
%                   Mueller matrix exceeds +/-1.
%   (28Jan2021) -Instead of finding the physical MM (PMM) of the normalized calibrated MM (NMM),
%                   now using the unnormalized calibrated MM to find PMM.
%   (01Feb2021) -Before being acted on by DRM, irradiances are now normalized by air
%                   calibration irradiances, rather than by itself.
%   (22Feb2021) -Discovered 'MuellerM' is incorrectly calculating Mueller
%                   matrix through mis-transposing elements; RESOLVED 22Feb2021
%               -Discovered 'PhysMM' function is incorrectly calculating
%                   physical Mueller matrix; RESOLVED 13Apr2021
%   (24Feb2021) -Added option to toggle whether to plot results.
%   (11Mar2021) -Added Mueller plot function 'PlotMueller'.
%               -Changed 'Intensity' to 'Irradiance', in name only.
%   (24Mar2021) -Included measurement mode 'reflection'.
%               -Gave 'IDIFF' its own plot whereas before it was plotted
%                   with the systemic calibration errors.
%   (29Mar2021) -Created function 'RefIDXcatalog' to provide the real and
%                   imaginary components of the complex refractive index of a material as a
%                   function of wavelength.
%   (06May2021) -Must now specify sample's 4x4 ideal Mueller matrix in
%                   sample's data folder. Not required for either transmission or 
%                   reflection calibration measurements (built-in).
%   (13May2021) -Must now create corresponding 'PARAMETERS' text file to go
%                   with measurement data; code reads text file for desired
%                   settings.
%   (08June2021) -Incorporated flexiblity into code to manage differing wavelength
%                   ranges between sample data and calibration data, as long as desired sample
%                   waveband is within calibration waveband.
%
% NOTES:
%   Functions specific to reflection mode have been removed from the file for the
%   github upload, since the vortex retarder measurements were performed
%   in transmission mode. -Ella James, 16 Feb.2024
%==========================================================================
tic
clear 
addpath('POLARfuncs')
close all
w = warning('off','all'); % suppress display of all warning messages (such as creating a folder when it...
                            ...already exists); error messages still appear though

%% FILES/FOLDER
%==========================================================================
% POLARIMETER DATA FILES
% Spectropolarimeter LabVIEW data files without '.lvm' extension
FL.Calibration.FLDR = 'Data\DiagnosticRunsTRANSMISSION\24Jan2023';
FL.Calibration.NAME = 'CAL1_24Jan2023'; 

FL.Sample.FLDR = 'Data\SamplesTRANSMISSION\VortexRetarderWPV10L532';
FL.Sample.NAME = 'WPV10L532_07Feb2023';

%% TOGGLES
%==========================================================================
% CALIBRATION TOGGLE
TOG.Cal = 0; % '1'-> use calibration file to calculate system's corrected 'DRM' and systemic errors;...
                    ...else apply corrected 'DRM' to sample data file

% COMPARE TOGGLE
TOG.Comp = 0; % '1'-> calculate and display difference of sample's calibrated/physical Mueller 
                    ...matrix from ideal, and display ideal MM on all MM plots; else don't

% PLOT TOGGLE
TOG.Plot = 1; % '1'-> plot results and save figures; else don't
          
%% MAIN FUNCTION
%==========================================================================
POLARmain(FL,TOG);
tf = toc;
