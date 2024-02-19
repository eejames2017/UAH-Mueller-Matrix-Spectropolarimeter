function [NMM,MM] = POLARmain(FL,TOG)

% Last Edited 07 June 2021
%==========================================================================
% DESCRIPTION:
%   Innards of POLAR.m script.
%
% CREATED BY: 
%   Prashant Raman, UAH Physics PhD Student, 2012
% 
% EDITED BY: 
%   Matthew Goforth, UAH Physics Master's Student, 17 Nov 2020
%==========================================================================

%% SETUP
%==========================================================================
% FILES / PATHS
if TOG.Cal == 1
    fldr = fullfile(FL.Calibration.FLDR,FL.Calibration.NAME); % folder with polarimeter data and MATLAB results
    DRMfldr = fldr; % folder with DRM     
    ext = '.lvm';
        flnm = sprintf('%s%s',FL.Calibration.NAME,ext); % tac on extension to filename
    DATAfl = fullfile(fldr,flnm); % calibration file name/path
else 
    fldr = fullfile(FL.Sample.FLDR,FL.Sample.NAME); % folder with polarimeter data and MATLAB results
    DRMfldr = fullfile(FL.Calibration.FLDR,FL.Calibration.NAME); % folder with DRM
    ext = '.lvm';
        flnm = sprintf('%s%s',FL.Sample.NAME,ext); % tac on extension to filename
    DATAfl = fullfile(fldr,flnm); % calibration file name/path
        CALfl = sprintf('%s%s',FL.Calibration.NAME,ext); % tac on extension to filename
    CALfl = fullfile(FL.Calibration.FLDR,FL.Calibration.NAME,CALfl);
end

% PARAMETERS/SETTINGS
PARAflnm = 'PARAMETERS.txt';
PARA = ReadParameters(PARAflnm,fldr);
MODE = PARA.MODE;
INCIDENTnum = str2double(PARA.INCIDENT(1:end-3));
PARA.INCIDENTnum = INCIDENTnum;
if strcmpi(MODE,'reflection') == 1
    TOG.Mode = 1;
else
    TOG.Mode = 0;
end

% SORT POLARIMETER DATA
% Separate raw polarimeter data into arrays of angles, wavelengths, and irradiances
if TOG.Cal == 1
    [ANG,WL,IRR,NumWL,~] = PreProcess(DATAfl);
    calIRR = IRR;
    calIDX = ismember(WL,WL);
else
    [~,calWL,calIRR,~,~] = PreProcess(CALfl);
    [ANG,WL,IRR,~,~] = PreProcess(DATAfl);
    
    % Adjust accordingly if calibration wavelengths are different from
    % sample wavelengths
    calIDX = ismember(calWL,WL);
    sampIDX = ismember(WL,calWL);
    
    WL = WL(sampIDX);
    IRR = IRR(sampIDX,:);
    calIRR = calIRR(calIDX,:);
    NumWL = length(WL);
end

%% SELECT/SPAN IDEAL MUELLER MATRIX
%==========================================================================
% Replace/select appropriate ideal Mueller matrix depending on measurement mode
% and if it is a calibration measurement; then span 'MMideal' across wavelength band
if TOG.Cal == 1
    if TOG.Mode == 1 % if performing calibration, change ideal MM based on measurement mode; otherwise keep ideal MM as it is
        MMideal = zeros(4,4,NumWL);
        [n,k] = RefIDXcatalog(PARA.REFIDX,WL);
        for ii = 1:NumWL
            MMideal(:,:,ii) = MMmirror(n(ii),k(ii),INCIDENTnum); % Mueller matrix of mirror used to calibrate operation in reflection mode; uses dispersive complex refractive index
        end
    else
        MMideal = repmat(eye(4),1,1,NumWL); % Mueller matrix of air is the identity matrix over most wavelength bands; used in transmission mode
    end
else
    load(fullfile(fldr,'MMideal.mat'),'MMideal'); 
    MMideal = repmat(MMideal,1,1,NumWL); % retain ideal MM specified in sample's raw data folder and span over wavelength band
end

% Normalize by 'm11'
NMMideal = zeros(4,4,NumWL);
for ii = 1:NumWL 
    NMMideal(:,:,ii) = MMideal(:,:,ii)/MMideal(1,1,ii);
end
PARA.MMideal = MMideal;
PARA.NMMideal = NMMideal;

%% PROCESS POLARIMETER DATA
%==========================================================================
if TOG.Cal == 1
    % CALIBRATION ERRORS
    % Determine system's calibration errors through least-squares fit
    [CalData] = Calibration(ANG,WL,IRR,fldr,TOG,PARA);

    % DATA REDUCTION MATRIX
    % Generate polarimeter's calibrated data reduction matrix 'DRM' using
    % calibration errors 'CalData'; also generate raw/uncalibrated data
    % reduction matrix 'RDRM' for comparison
    [~,~] = DRmatrix(ANG,WL,CalData,fldr);
end

% MUELLER MATRIX
% Calculate & plot sample's normalized Mueller matrix 'NMM' using 'DRM';
% also plot uncalibrated normalized Mueller matrix 'RNMM'. All but first 
% Mueller element normalized as mij/m11; m11 unormalized
[NMM,MM] = MuellerM(ANG,WL,IRR,calIRR,DRMfldr,fldr,TOG,PARA,calIDX);

% PHYSICAL MUELLER MATRIX
% Calculate MM's nearest normalized physical Mueller matrix 'NPMM'; plot 'NPMM'
[NPMM,~] = PhysMM(WL,MM,fldr,TOG,PARA);

% DECOMPOSE MUELLER MATRIX
% Deconstruct sample's normalized Mueller matrix (MM/m00) into constituent polarization element matrices
[dCOMP] = dcompMM(WL,NPMM,fldr);

% POLARIZATION PROPERTIES
% Calculate (linear, circular) retardance, diattenuation,..., and plot results
[~] = PolarProp(WL,dCOMP,fldr,TOG,PARA);
end
