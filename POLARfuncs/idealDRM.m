function [W,DRM] = idealDRM(Q)

% Last Edited 13 Apr 2021
%==========================================================================
% DESCRIPTION:
%   16xQ polarimetric data reduction matrix for dual rotating retarder
%   polarimeter for 'Q' measurements. The retarders are rotated in angular
%   incremenets of 5*dth:dth (R2:R1), where 'dth' is the angular increment
%   of R1.
% 
% REFERENCE:
%   -Chpt. 22 'Hanbook of Optics', vol.II, 2nd ed., pg.22.21
%   -Chenault, Pezzaniti, Chipman, 1992.
%
% CREATED BY:
%   Matt Goforth, 13 Apr 2021
%==========================================================================

%% PARAMETERS
%==========================================================================
d1 = 90; % R1 retardance (deg)
d2 = 90; % R2 retardance (deg)
dth = 360/(Q-1); % R1 angular increment (deg)

%% CONSTRUCT FLATTENED POLARIMETRIC MEASUREMENT MATRIX
%==========================================================================
% Precursor to the data reduction matrix.

W = zeros(Q,16);
for ii = 1:Q
    q = ii-1;
    
    % ELEMENTS OF qth FLATTENED POLARIMETRIC MEASUREMENT MATRIX
    wq1 = 1;
    wq2 = cosd(d1/2)^2 + sind(d1/2)^2*cosd(4*q*dth);
    wq3 = sind(d1/2)^2*sind(4*q*dth);
    wq4 = sind(d1)*sind(2*q*dth);

    wq5 = cosd(d2/2)^2 + sind(d2/2)^2*cosd(20*q*dth);
    wq6 = cosd(d1/2)^2*cosd(d2/2)^2 + sind(d1/2)^2*cosd(d2/2)^2*cosd(4*q*dth) +...
          cosd(d1/2)^2*sind(d2/2)^2*cosd(20*q*dth) +...
          1/2*sind(d1/2)^2*sind(d2/2)^2*(cosd(16*q*dth) + cosd(24*q*dth));
    wq7 = sind(d1/2)^2*cosd(d2/2)^2*sind(4*q*dth) +...
          1/2*sind(d1/2)^2*sind(d2/2)^2*(-sind(16*q*dth) + sind(24*q*dth));
    wq8 = sind(d1)*cosd(d2/2)^2*sind(2*q*dth) +...
          1/2*sind(d1)*sind(d2/2)^2*(-sind(18*q*dth) + sind(22*q*dth));

    wq9 = sind(d2/2)^2*sind(20*q*dth);
    wq10 = cosd(d1/2)^2*sind(d2/2)^2*sind(20*q*dth) +...
           1/2*sind(d1/2)^2*sind(d2/2)^2*(sind(16*q*dth) + sind(24*q*dth));
    wq11 = 1/2*sind(d1/2)^2*sind(d2/2)^2*(cosd(16*q*dth) - cosd(24*q*dth));
    wq12 = 1/2*sind(d1)*sind(d2/2)^2*(cosd(18*q*dth) - cosd(22*q*dth));
    
    wq13 = -sind(d2)*sind(10*q*dth);
    wq14 = -cosd(d1/2)^2*sind(d2)*sind(10*q*dth) -...
            1/2*sind(d1/2)^2*sind(d2)*(sind(6*q*dth) + sind(14*q*dth));
    wq15 = -1/2*sind(d1/2)^2*sind(d2)*(cosd(6*q*dth) - cosd(14*q*dth));
    wq16 = -1/2*sind(d1)*sind(d2)*(cosd(8*q*dth) - cosd(12*q*dth));

    % FLATTENED POLARIMETRIC MEASUREMENT MATRIX
    W(ii,:) = 0.25*[wq1 wq2 wq3 wq4 wq5 wq6 wq7 wq8 wq9 wq10 wq11 wq12 wq13 wq14 wq15 wq16];
end

W = W(1:end-1,:);

%% DATA REDUCTION MATRIX
%==========================================================================
DRM = pinv(W);

end

