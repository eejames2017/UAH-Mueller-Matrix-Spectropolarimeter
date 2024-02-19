function MM = MMret(phi,theta)

% Last Edited 25 Mar 2021
%==========================================================================
% DESCRIPTION:
%   Mueller matrix for ideal retarder having retardance 'phi' and
%   fast axis 'theta', both in degrees. 'theta = 0deg' corresponds to the
%   fast axis being horizontally aligned with and in the direction of the
%   positive x-axis.
% 
% REFERENCE:
%   Definition for MM of retarder taken from Absorption and Scattering of...
%   Light by Small Particles, Bohren & Huffman, pg 55
%
% CREATED BY:
%   Matt Goforth, 09 Sep 2020
%==========================================================================

%% MUELLER ELEMENTS
%==========================================================================
m22 = (cosd(2*theta)^2)+(cosd(phi)*(sind(2*theta)^2));
m23 = (1-cosd(phi))*cosd(2*theta)*sind(2*theta);
m24 = -sind(2*theta)*sind(phi);

m32 = (1-cosd(phi))*cosd(2*theta)*sind(2*theta);
m33 = (sind(2*theta)^2)+(cosd(phi)*(cosd(2*theta)^2));
m34 = cosd(2*theta)*sind(phi);

m42 = sind(2*theta)*sind(phi);
m43 = -cosd(2*theta)*sind(phi);
m44 = cosd(phi);

%% MUELLER MATRIX
%==========================================================================
MM = [1   0   0   0;
      0 m22 m23 m24;
      0 m32 m33 m34;
      0 m42 m43 m44];
end

