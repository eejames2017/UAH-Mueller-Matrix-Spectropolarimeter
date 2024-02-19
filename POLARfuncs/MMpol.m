function MM = MMpol(theta)

% Last Edited 25 Mar 2021
%==========================================================================
% DESCRIPTION:
%   Mueller matrix for linear polarizer rotated 'theta' degrees.
% 
% REFERENCE:
%  ...
%
% EDITED BY:
%   Matt Goforth, 09 Sept 2020
%==========================================================================

%% PARAMETERS
%==========================================================================
% Amplitude attenuation coefficient
p = 1; % p = 1 for an ideal linear polarizer

% Quality of polarizer (related to extinction ratio)
gamma = 0; % Ranges from 0 <= gamma <= 90deg, where...
           ...gamma = 0 deg --> ideal horizontal polarizer
           ...gamma = 90 deg --> ideal vertical polarizer

%% MUELLER ELEMENTS
%==========================================================================
m12 = cosd(2*gamma)*cosd(2*theta);
m13 = cosd(2*gamma)*sind(2*theta);

m21 = cosd(2*gamma)*cosd(2*theta);
m22 = cosd(2*theta)^2 + sind(2*gamma)*(sind(2*theta)^2);
m23 = (1-sind(2*gamma))*cosd(2*theta)*sind(2*theta);

m31 = cosd(2*gamma)*sind(2*theta);
m32 = (1-sind(2*gamma))*cosd(2*theta)*sind(2*theta);
m33 = sind(2*theta)^2 + sind(2*gamma)*(cosd(2*theta)^2);

m44 = sind(2*gamma);

%% MUELLER MATRIX
%==========================================================================
MM = (p^2)*0.5*[  1 m12 m13   0; 
                m21 m22 m23   0; 
                m31 m32 m33   0; 
                  0   0   0 m44];
end

