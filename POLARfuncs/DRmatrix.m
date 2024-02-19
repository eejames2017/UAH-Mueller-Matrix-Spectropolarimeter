function [RDRM,DRM] = DRmatrix(Angles,Wavelengths,CalData,fldr)

% Last Edited 07 June 2021
%==========================================================================
% DESCRIPTION:
%   Calculates data reduction matrix for measurement, accounting for R1 and
%   R2 orientation errors and R1, R2, P2 orientation errors
% 
% REFERENCE:
% P. Raman, “Spectropolarimetric characterization of light scattering materials,” 
% Ph.D. thesis, University of Alabamain Huntsville (2012): Chapter 3
%
% EDITED BY:
%   Matt Goforth, 02 Nov 2020
%==========================================================================

%% SETUP
%==========================================================================
% READ IN CALIBRATION DATA
e1 = CalData(:,1); % retardance error of R1 relative to ideal quarter-wave retardance (deg)
e2 = CalData(:,2); % retardance error of R2 relative to ideal quarter-wave retardance (deg)
e3 = CalData(:,3); % orientation error of R1 relative to fast-axis of P1 (deg)
e4 = CalData(:,4); % orientation error of R2 '...' (deg)
e5 = CalData(:,5); % orientation error of P2 '...' (deg)

% MISCELLANEOUS   
NumWL = length(Wavelengths); % number of wavelengths scanned over by polarimeter
NumANG = length(Angles); % number of angles assumed by R1 & R2 during each measurement sequence

tol = 10^(-10); % tolerance for an element of W matrix to be considered essentially 0 (see below)

%% CALCULATIONS
%==========================================================================
% PREALLOCATION
dR1 = zeros(NumWL,NumANG);
dR2 = zeros(NumWL,NumANG);
thR1 = zeros(NumWL,NumANG);
thR2 = zeros(NumWL,NumANG);
thP2 = zeros(NumWL,NumANG);

RW = zeros(NumANG,16,NumWL);
RDRM = zeros(16,NumANG,NumWL); 
W = zeros(NumANG,16,NumWL);
DRM = zeros(16,NumANG,NumWL); % since DRM is the pseudo inverse of W, first two dimensions are swapped

for ii = 1:NumWL
    dR1(ii,:) = 90+e1(ii); % corrected retardance of R1
    dR2(ii,:) = 90+e2(ii); % corrected retardance of R2
    thR1(ii,:) = Angles+e3(ii); % corrected orientation of R1
    thR2(ii,:) = (5*Angles)+e4(ii); % corrected orientation of R2
    thP2(ii,:) = e5(ii); % corrected orientation of P2
end

% UNCALIBRATED DATA REDUCTION MATRIX
for ii = 1:NumWL
    for jj = 1:NumANG
        Sgenerator = MMret(90,Angles(jj))*[1 1 0 0]'; % Stokes vector produced from polarization generator optical elements
        MManalyzer = MMpol(0)*MMret(90,5*Angles(jj)); % Mueller matrix of the polarization analyzer optical elements
        Avect = MManalyzer(1,:)'; % analyzer vector
        RW(jj,:,ii) = reshape((Avect*Sgenerator')',1,16); % polarimetric measurement matrix for uncalibrated DRM     
    end
    RW(abs(RW) < tol) = 0; % replace elements of the polarimetric measurement matrix with 0 if < tolerance
    RDRM(:,:,ii) = pinv(RW(:,:,ii), tol); % raw data reduction matrix
    
end

% CALIBRATED DATA REDUCTION MATRIX
for ii = 1:NumWL
    for jj = 1:NumANG
        Sgenerator = MMret(dR1(ii,jj),thR1(ii,jj))*[1 1 0 0]'; % Stokes vector produced from polarization generator optical elements
        MManalyzer = MMpol(thP2(ii,jj))*MMret(dR2(ii,jj),thR2(ii,jj)); % Mueller matrix of the polarization analyzer optical elements
        Avect = MManalyzer(1,:)'; % analyzer vector
        W(jj,:,ii) = reshape((Avect*Sgenerator')',1,16); % polarimetric measurement matrix       
    end
    W(abs(W) < tol) = 0; % replace elements of the polarimetric measurement matrix with 0 if < tolerance
    DRM(:,:,ii) = pinv(W(:,:,ii), tol); % data reduction matrix   
end

%% SAVE
%==========================================================================
flpthRDRM = fullfile(fldr,'RDRM');
save(flpthRDRM,'RDRM');

flpthDRM = fullfile(fldr,'DRM');
save(flpthDRM,'DRM');

end
