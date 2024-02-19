function [NPMM,PMM] = PhysMM(Wavelengths,MM,fldr,TOG,PARA)

% Last Edited 07 June 2021
%==========================================================================
% DESCRIPTION:
%   Improve the sample's calibrated Mueller matrix by finding its closest
%   physically realizable Mueller matrix through a process that removes 
%   negative eigenvalues/vectors.
% 
% REFERENCE:
%   S.Y. Lu, "An Interpretation of Polarization Matrices," Ph.D. Thesis,
%   University of Alabama in Huntsville (1995).
%
%   P. Raman, "Spectropolarimetric characterization of light scattering materials,” 
%   Ph.D. thesis, University of Alabama in Huntsville (2012): Chapter 3
%
% EDITED BY:
%   Matt Goforth, 06 Nov 2020
%==========================================================================

%% SETUP
%==========================================================================
% MISCELLANEOUS
NumWL = length(Wavelengths);

% N-MATRIX
% Define a complex Hermitian matrix 'Nmatrix' for each 4x4 real matrix 'MM',
% called its N-matrix. (pg 81 Shih-Yau Lu)
Nmatrix(1,1,:) = 0.5*(MM(1,1,:)+ MM(2,2,:)+    MM(1,2,:)+    MM(2,1,:));
Nmatrix(1,2,:) = 0.5*(MM(1,3,:)+ MM(2,3,:)+ 1i*MM(1,4,:)+ 1i*MM(2,4,:));
Nmatrix(1,3,:) = 0.5*(MM(3,1,:)+ MM(3,2,:)- 1i*MM(4,1,:)- 1i*MM(4,2,:));
Nmatrix(1,4,:) = 0.5*(MM(3,3,:)+ MM(4,4,:)+ 1i*MM(3,4,:)- 1i*MM(4,3,:));
Nmatrix(2,1,:) = 0.5*(MM(1,3,:)+ MM(2,3,:)- 1i*MM(1,4,:)- 1i*MM(2,4,:));
Nmatrix(2,2,:) = 0.5*(MM(1,1,:)- MM(2,2,:)-    MM(1,2,:)+    MM(2,1,:));
Nmatrix(2,3,:) = 0.5*(MM(3,3,:)- MM(4,4,:)- 1i*MM(3,4,:)- 1i*MM(4,3,:));
Nmatrix(2,4,:) = 0.5*(MM(3,1,:)- MM(3,2,:)- 1i*MM(4,1,:)+ 1i*MM(4,2,:));
Nmatrix(3,1,:) = 0.5*(MM(3,1,:)+ MM(3,2,:)+ 1i*MM(4,1,:)+ 1i*MM(4,2,:));
Nmatrix(3,2,:) = 0.5*(MM(3,3,:)- MM(4,4,:)+ 1i*MM(3,4,:)+ 1i*MM(4,3,:));
Nmatrix(3,3,:) = 0.5*(MM(1,1,:)- MM(2,2,:)+    MM(1,2,:)-    MM(2,1,:));
Nmatrix(3,4,:) = 0.5*(MM(1,3,:)- MM(2,3,:)+ 1i*MM(1,4,:)- 1i*MM(2,4,:));
Nmatrix(4,1,:) = 0.5*(MM(3,3,:)+ MM(4,4,:)- 1i*MM(3,4,:)+ 1i*MM(4,3,:));
Nmatrix(4,2,:) = 0.5*(MM(3,1,:)- MM(3,2,:)+ 1i*MM(4,1,:)- 1i*MM(4,2,:));
Nmatrix(4,3,:) = 0.5*(MM(1,3,:)- MM(2,3,:)- 1i*MM(1,4,:)+ 1i*MM(2,4,:));
Nmatrix(4,4,:) = 0.5*(MM(1,1,:)+ MM(2,2,:)-    MM(1,2,:)-    MM(2,1,:));
 
% TRANSFORMATION MATRIX
% transformation matrix from flattened Mueller matrix to flattened N-matrix
C = [1  1  0   0  1  1  0   0  0  0  0   0    0    0   0  0; 
     0  0  1  1i  0  0  1  1i  0  0  0   0    0    0   0  0; 
     0  0  0   0  0  0  0   0  1  1  0   0  -1i  -1i   0  0; 
     0  0  0   0  0  0  0   0  0  0  1  1i    0    0 -1i  1; 
     0  0  1 -1i  0  0  1 -1i  0  0  0   0    0    0   0  0; 
     1 -1  0   0  1 -1  0   0  0  0  0   0    0    0   0  0; 
     0  0  0   0  0  0  0   0  0  0  1 -1i    0    0 -1i -1; 
     0  0  0   0  0  0  0   0  1 -1  0   0  -1i   1i   0  0; 
     0  0  0   0  0  0  0   0  1  1  0   0   1i   1i   0  0; 
     0  0  0   0  0  0  0   0  0  0  1  1i    0    0  1i -1; 
     1  1  0   0 -1 -1  0   0  0  0  0   0    0    0   0  0; 
     0  0  1  1i  0  0 -1 -1i  0  0  0   0    0    0   0  0; 
     0  0  0   0  0  0  0   0  0  0  1 -1i    0    0  1i  1; 
     0  0  0   0  0  0  0   0  1 -1  0   0   1i  -1i   0  0; 
     0  0  1 -1i  0  0 -1  1i  0  0  0   0    0    0   0  0; 
     1 -1  0   0 -1  1  0   0  0  0  0   0    0    0   0  0];
 
% PREALLOCATE
estNmatrix = zeros(size(Nmatrix));
SumNegEIGsq = zeros(NumWL);
SumPosEIGsq = zeros(NumWL);
estERR = zeros(1,NumWL);
rmsERR = zeros(1,NumWL);

%% CALCULATIONS
%==========================================================================
% Process for finding calibrated Mueller matrix's nearest physically realizable
% Muller matrix taken from Chpt 5 of dissertation 'An Interpretation of Polarization Matrices' by Shih-Yau Lu.

% Criterion for physical Mueller matrices: A 4x4 real matrix is a physical
% Mueller matrix if and only if its N-matrix is nonnegative definite. (pg 90 Shih-Yau Lu)

% NEAREST PHYSICAL N-MATRIX & ESTIMATED ERROR
% Nearest physicsal N-matrix (estN-matrix) found by removing negative eigenvalues from
% N-matrix; estimated error 'estERR' quantifies how well estPhysN-matrix approximates
% N-matrix. (pg 91,94 Shih-Yau Lu)
for jj = 1:NumWL
    [eigVEC,eigVAL] = eig(Nmatrix(:,:,jj),'vector'); % gives eigenvalues in column vector 'eigVAL', and gives matrix 'eigVEC' whose...
                                                    ...columns are the corresponding right eigenvectors
    for kk = 1:4
        if eigVAL(kk) >= 0
            estNmatrix(:,:,jj) = estNmatrix(:,:,jj)+eigVAL(kk)*(eigVEC(:,kk)*eigVEC(:,kk)'); % Eq 5.30 Shih-Yau Lu
            SumPosEIGsq(jj) = SumPosEIGsq(jj)+eigVAL(kk)^2;
        else
            SumNegEIGsq(jj) = SumNegEIGsq(jj)+eigVAL(kk)^2;
        end
    end
    estERR(jj) = 100*sqrt(SumNegEIGsq(jj)/(SumNegEIGsq(jj)+SumPosEIGsq(jj)));
end

% PHYSICAL MUELLER MATRIX
PMM = zeros(4,4,NumWL); % preallocate
NPMM = zeros(4,4,NumWL); % preallocate
for jj = 1:NumWL
    X = reshape(estNmatrix(:,:,jj).',[16,1]); % estNmatrix().' is the nonconjugate transpose of estNmatrix()
    M = C\X; % solves linear equation 'C*M = X' for 'M'
    PMM(:,:,jj) = reshape(M,4,4).';
    PMM(:,:,jj) = 2*real(PMM(:,:,jj)); % take real components of 'PMM'
    NPMM(:,:,jj) = PMM(:,:,jj)/PMM(1,1,jj); % normalize 'PMM'
end

% SANITY CHECK
% Terminate if... 
if abs(NPMM) > 1
    error('An Element in Normalized Physical Mueller Matrix Exceeds %c1!', char(177)); % ...a Mueller element exceeds +/-1
end

if TOG.Comp == 1
    % DIFFERENCE OF NORMALIZED PHYSICAL MUELLER MATRIX FROM NORMALIZED IDEAL
    NPMMdiff = NPMM - PARA.NMMideal; 
    NPMMdiff(1,1,:) = PMM(1,1,:) - PARA.MMideal(1,1,:);

    % RMS ERROR
    % Frobenius norm of difference between nearest physical MM measured and ideal MM
    for jj = 1:NumWL
        rmsERR(jj) = 100*norm(NPMMdiff(:,:,jj),'fro')/4;
    end 
    
end

%% PLOT 
%==========================================================================
if TOG.Plot == 1
    % SETUP
    fWL = Wavelengths(1);
    eWL = Wavelengths(end);
    key = reshape(1:16,[4 4])'; % place Mueller element in appropriate subplot
                                % position given indices 'ii' & 'jj'
    
    % PLOT NORMALIZED PHYSICAL MUELLER MATRIX
    fig1 = figure();
    TITLE1 = 'Normalized Physical Mueller Matrix';
    PlotMueller(Wavelengths,PMM,NPMM,TITLE1,fig1,PARA,TOG.Comp,0)
    
    % PLOT ESTIMATED ERROR
    fig2 = figure();
    plot(Wavelengths,estERR,'-r','LineWidth',1)
    xlabel('\lambda (nm)','FontSize',12)
    title(sprintf('%%RMSerror of Measured PMM from Measured MM,  %s',PARA.SAMPLE),'FontSize',20,'FontWeight','Normal')
    xlim([fWL eWL])
    grid on
    grid minor
    ax1 = gca;
    
    hold on
    pa = plot(nan,nan,'ko',nan,nan,'ko',nan,nan,'ko',nan,nan,'ko',nan,nan,'ko',nan,nan,'ko'); % create 2nd legend for listing lab setup parameters
    hold off

    ax1POSITION = ax1.Position; % position of first axes
    ax2 = axes('position',ax1POSITION,'visible','off'); % create 2nd invisible axes for 2nd legend
    LEG2 = legend(ax2,pa,...
            sprintf('%s',PARA.MODE),...
            sprintf('{\\Delta\\lambda}_{Resolution} = %s',PARA.RESOLUTION),...
            sprintf('{\\Delta\\lambda}_{Step} = %s',PARA.STEP),...
            sprintf('{\\lambda}_{o} = %s',PARA.LAMBDAo),... 
            sprintf('{\\lambda}_{f} = %s',PARA.LAMBDAf),...
            sprintf('{\\theta}_{incident} = %s',PARA.INCIDENT)); 
    LEG2.Title.String = sprintf('%s',PARA.DATE);
    LEG2.Position = [0.00184688758419097 0.646565540252065 0.0911458314086001 0.260963082137033];
    LEG2.FontSize = 11;  
    
    set(fig2,'Position', get(0,'Screensize')); % maximizes figure window; looks better when saving figure

    if TOG.Comp == 1
        % PLOT DIFFERENCE OF NPMM FROM NORMALIZED IDEAL MM
        fig3 = figure();
        for ii = 1:4
            for jj = 1:4
                subplot(4,4,key(ii,jj))
                plot(Wavelengths,squeeze(NPMMdiff(ii,jj,:)),'r','LineWidth',1)
                if ii == 1 && jj == 1
                    ylabel(sprintf('m_{%u%u} - m^{Ideal}_{%u%u}',ii,jj,ii,jj),'FontSize',12)
                else
                    ylabel(sprintf('m_{%u%u}/m_{11} - m^{Ideal}_{%u%u}/m^{Ideal}_{11}',ii,jj,ii,jj),'FontSize',12)
                end
                xlabel('\lambda (nm)','FontSize',12)
                title(sprintf('m_{%u%u}',ii,jj),'FontSize',12)
                xlim([fWL eWL])
                grid on
                box on   
            end
        end
        sgtitle(sprintf('Difference between Measured NPMM and Ideal NMM,  %s',PARA.SAMPLE),'FontSize',20)
        ax1 = gca;
    
        hold on
        pa = plot(nan,nan,'ko',nan,nan,'ko',nan,nan,'ko',nan,nan,'ko',nan,nan,'ko',nan,nan,'ko'); % create 2nd legend for listing lab setup parameters
        hold off

        ax1POSITION = ax1.Position; % position of first axes
        ax2 = axes('position',ax1POSITION,'visible','off'); % create 2nd invisible axes for 2nd legend
        LEG2 = legend(ax2,pa,...
                sprintf('%s',PARA.MODE),...
                sprintf('{\\Delta\\lambda}_{Resolution} = %s',PARA.RESOLUTION),...
                sprintf('{\\Delta\\lambda}_{Step} = %s',PARA.STEP),...
                sprintf('{\\lambda}_{o} = %s',PARA.LAMBDAo),... 
                sprintf('{\\lambda}_{f} = %s',PARA.LAMBDAf),...
                sprintf('{\\theta}_{incident} = %s',PARA.INCIDENT)); 
        LEG2.Title.String = sprintf('%s',PARA.DATE);
        LEG2.Position = [0.00184688758419097 0.646565540252065 0.0911458314086001 0.260963082137033];
        LEG2.FontSize = 11;  
        
        set(fig3,'Position', get(0,'Screensize')); % maximizes figure window; looks better when saving figure

        % PLOT RMS ERROR
        fig4 = figure();
        plot(Wavelengths,rmsERR,'-r','LineWidth',1);
        xlabel('\lambda (nm)','FontSize',12);
        title(sprintf('%%RMSerror of Measured NPMM from Ideal NMM,  %s',PARA.SAMPLE),'FontSize',20,'FontWeight','Normal');
        xlim([fWL eWL])
        grid on
        grid minor
        ax1 = gca;
    
        hold on
        pa = plot(nan,nan,'ko',nan,nan,'ko',nan,nan,'ko',nan,nan,'ko',nan,nan,'ko',nan,nan,'ko'); % create 2nd legend for listing lab setup parameters
        hold off

        ax1POSITION = ax1.Position; % position of first axes
        ax2 = axes('position',ax1POSITION,'visible','off'); % create 2nd invisible axes for 2nd legend
        LEG2 = legend(ax2,pa,...
                sprintf('%s',PARA.MODE),...
                sprintf('{\\Delta\\lambda}_{Resolution} = %s',PARA.RESOLUTION),...
                sprintf('{\\Delta\\lambda}_{Step} = %s',PARA.STEP),...
                sprintf('{\\lambda}_{o} = %s',PARA.LAMBDAo),... 
                sprintf('{\\lambda}_{f} = %s',PARA.LAMBDAf),...
                sprintf('{\\theta}_{incident} = %s',PARA.INCIDENT)); 
        LEG2.Title.String = sprintf('%s',PARA.DATE);
        LEG2.Position = [0.00184688758419097 0.646565540252065 0.0911458314086001 0.260963082137033];
        LEG2.FontSize = 11;  
        
        set(fig4,'Position', get(0,'Screensize')); % maximizes figure window; looks better when saving figure
    end
end

%% SAVE PLOTS
%==========================================================================
flpthEstERR = fullfile(fldr,'estERR'); 
save(flpthEstERR,'estERR'); % save EstERR.mat data file

flpthNPMM = fullfile(fldr,'NPMM'); 
save(flpthNPMM,'NPMM'); % save NPMM.mat data file

flpthPMM = fullfile(fldr,'PMM'); 
save(flpthPMM,'PMM'); % save PMM.mat data file

if TOG.Plot == 1
    saveas(fig2,flpthEstERR,'jpg'); % save EstErr plot as .jpg
    saveas(fig1,flpthNPMM,'jpg'); % save NPMM plot as .jpg
end

if TOG.Comp == 1
    flpthNPMMdiff = fullfile(fldr,'NPMMdiff'); 
    save(flpthNPMMdiff,'NPMMdiff'); % save NPMMdiff.mat data file
    
    flpthRMSerr = fullfile(fldr,'rmsERR'); 
    save(flpthEstERR,'rmsERR'); % save rmsERR.mat data file
    if TOG.Plot == 1
        saveas(fig3,flpthNPMMdiff,'jpg'); % save NPMMdiff plot as .jpg
        saveas(fig4,flpthRMSerr,'jpg'); % save rmsERR plot as .jpg
    end
end
end
