function normInt = ModelFunctionTRANS(error,theta)

% Last Edited 25 Mar 2021
%==========================================================================
% DESCRIPTION:
%   Calculates normalized irradiance modulation, correcting for R1 and R2
%   retardance errors and R1, R2, P2 orientation errors. The irradiances
%   here are the "optimal" values used the least-square solution for the
%   calibration errors. This function is called iteratively until the
%   normalized irradiance from least-squares solution for the 5 errors 
%   matches what is calculated here.
% 
% 
% REFERENCE:
%   P. Raman, "Spectropolarimetric characterization of light scattering materials,” 
%   Ph.D. thesis, University of Alabama in Huntsville (2012) : Chapter 3
%
% EDITED BY:
%   Matt Goforth, 30 Oct 2020
%==========================================================================

%% ERRORS
%==========================================================================
e1 = error(1); % retardance error of R1 relative to ideal quarter-wave plate (deg)
e2 = error(2); % retardance error of R2 relative to ideal quarter-wave plate (deg)
e3 = error(3); % orientation error of R1 relative to fast-axis of P1 (deg)
e4 = error(4); % orientation error of R2 '...' (deg)
e5 = error(5); % orientation error of P2 '...' (deg)

%% CORRECTED PARAMETERS
%==========================================================================
dR1 = 90+e1; % corrected retardance of R1 (deg) 
dR2 = 90+e2; % corrected retardance of R2 (deg)
thR1 = theta+e3; % corrected orientation of R1 (deg)
thR2 = (5*theta)+e4; % corrected orientation of R2 (deg)
thP2 = e5; % corrected orientation of P2 (deg)

%% CALCULATED NORMALIZED IRRADIANCE
%==========================================================================
Irr = zeros(1,length(theta)); % preallocate
for ii = 1:length(theta)
    StokesOut = MMpol(thP2)*MMret(dR2,thR2(ii))*MMair()*MMret(dR1,thR1(ii))*[1 1 0 0]';
    Irr(ii) = StokesOut(1);
end

normInt = Irr/max(Irr); % normalize irradiances

end