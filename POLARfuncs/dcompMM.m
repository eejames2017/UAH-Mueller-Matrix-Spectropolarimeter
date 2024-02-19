function [dCOMP] = dcompMM(Wavelengths,NMM,fldr)

% Last Edited 07 June 2021
%==========================================================================
% DESCRIPTION:
%   Decompose normalized Mueller matrix 'NMM' into its constituent matrices each
%   representing a different polarizing element (retarder, diattenuater,
%   depolarizer).
% 
% REFERENCE:
%   S. Y. Lu and R. A. Chipman, “Interpretation of Mueller matrices based on polar decomposition,” J. Opt. Soc. Am. A 13, 1106–1113 (1996).
%
% EDITED BY:
%   Matt Goforth, 06 Nov 2020
%==========================================================================

%% SETUP
%==========================================================================
% MISCELLANEOUS
NumWL = length(Wavelengths);

% DIATTENUATION VECTOR AND POLARIZANCE VECTOR
Dvect = reshape(NMM(1,2:4,:),3,NumWL); % diattenuation vector
Pvect = reshape(NMM(2:4,1,:),3,NumWL); % polarizance vector

%% CONSTITUENT POLARIZING ELEMENT MATRICES
%==========================================================================
% CONSTITUENT DIATTENUATER MUELLER MATRIX
Tu = 1; % transmittance for incident unpolarized light; not sure what this should be, M.G. 02 June 2021
D = sqrt(dot(Dvect,Dvect)); % diattenuation

mmD = zeros(3,3,NumWL); % sub diattenuater Mueller matrix
for ii = 1:NumWL
    if D(ii) < 10^-6
        mmD(:,:,ii) = eye(3);
    else
        unitDvect = Dvect(:,ii)/D(ii);
        mmD(:,:,ii) = (sqrt(1-D(ii)^2)*eye(3))+(1-sqrt(1-D(ii)^2))*(unitDvect*unitDvect.');
    end
end

MMD = ones(4,4,NumWL); 
MMD(1,2:4,:) = Dvect;
MMD(2:4,1,:) = Dvect;
MMD(2:4,2:4,:) = mmD;
MMD = Tu*MMD; % diattenuater Mueller matrix

% CONSTITUENT DEPOLARIZATION MUELLER MATRIX
NoDiat = zeros(4,4,NumWL); % NoDiat is the M'matrix in Lu & Chipman paper (1996).
for ii = 1:NumWL
    NoDiat(:,:,ii) = NMM(:,:,ii)/MMD(:,:,ii);
end

DePvect = reshape(NoDiat(2:4,1,:),3,NumWL);
NoDiatSubMAT = NoDiat(2:4,2:4,:);
DetNoDiat = zeros(1,NumWL);
SignDet = zeros(1,NumWL);
eVAL = zeros(3,NumWL);
SqNoDiatSubMAT = zeros(3,3,NumWL);
Depolsubmat = zeros(3,3,NumWL);

for ii = 1:NumWL
    SqNoDiatSubMAT(:,:,ii) = NoDiatSubMAT(:,:,ii)*NoDiatSubMAT(:,:,ii)';
    DetNoDiat(ii) = det(NoDiatSubMAT(:,:,ii));
    if DetNoDiat(ii) > 0
        SignDet(ii) = 1;
    elseif DetNoDiat(ii) < 0
        SignDet(ii) = -1;
    else
        SignDet(ii) = 0;
    end
    eVAL(:,ii) = eig(SqNoDiatSubMAT(:,:,ii)');
    
    SQRT1 = sqrt(eVAL(1,ii)*eVAL(2,ii))+sqrt(eVAL(2,ii)*eVAL(3,ii))+sqrt(eVAL(3,ii)*eVAL(1,ii));
    SQRT2 = sqrt(eVAL(1,ii))+sqrt(eVAL(2,ii))+sqrt(eVAL(3,ii));
    SQRT3 = sqrt(eVAL(1,ii)*eVAL(2,ii)*eVAL(3,ii));
    MAT1 = SqNoDiatSubMAT(:,:,ii)+SQRT1*eye(3);
    MAT2 = SQRT2*SqNoDiatSubMAT(:,:,ii)+SQRT3*eye(3);
    
    Depolsubmat(:,:,ii) = SignDet(ii)*(MAT1\MAT2); 
end


depolMM = zeros(4,4,NumWL); % depolarizer Mueller matrix
temp = zeros(4,NumWL);

temp(1,:) = 1;
depolMM(1,:,:) = temp;
depolMM(2:4,1,:) = DePvect;
depolMM(2:4,2:4,:) = Depolsubmat;


% CONSTITUENT RETARDER MUELLER MATRIX
MMR = zeros(4,4,NumWL);
for ii = 1:NumWL
    MMR(:,:,ii) = depolMM(:,:,ii)\NoDiat(:,:,ii);
end

% CALCULATE RETARDANCE
R = zeros(1,NumWL); % retardance
for ii = 1:NumWL
    R(ii) = acosd(trace(MMR(:,:,ii))/2-1);
    if ii == 281 || ii == 282
        disp('R')
        disp(ii)
    disp(R(ii))
    end
end

% CALCULATE HORIZONTAL, 45-OBLIQUE, CIRCULAR RETARDANCE
a = zeros(3,NumWL); % retardance fast-axis vector (a1 -> RH, a2 -> R45, a3 -> Rcircular)
Rvect = zeros(3,NumWL); % retardance vector
for ii = 1:NumWL
    for jj = 1:3
        for kk = 1:3
            for mm = 1:3
                a(jj,ii) = a(jj,ii)+(((jj-kk)*(kk-mm)*(mm-jj)/2)*MMR(kk+1,mm+1,ii)/(2*sind(R(ii))));
            end
        end
    end
    Rvect(1,ii) = R(ii)*a(1,ii);
    Rvect(2,ii) = R(ii)*a(2,ii);
    Rvect(3,ii) = R(ii)*a(3,ii);
end


%% SAVE
%==========================================================================
% SORT INTO STRUCTURE
dCOMP.Dvect = Dvect;
dCOMP.Pvect = Pvect;
dCOMP.Rvect = Rvect;
dCOMP.depolMM = depolMM;

% SAVE TO FOLDER
save(fullfile(fldr,'dCOMP'),'dCOMP'); % save dCOMP structure
end
