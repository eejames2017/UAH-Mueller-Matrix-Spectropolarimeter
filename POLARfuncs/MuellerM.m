function [NMM,MM] = MuellerM(Angles,Wavelengths,Irradiances,CALIrradiances,DRMfldr,fldr,TOG,PARA,calIDX)

% Last Edited 07 June 2021
%==========================================================================
% DESCRIPTION:
%   Calculates the Mueller matrix of the sample
% 
% REFERENCE:
%   ...
%
% EDITED BY:
%   Matt Goforth, 02 Nov 2020
%==========================================================================

%% SETUP
%==========================================================================
% LOAD-IN DATA REDUCTION MATRICES
load(fullfile(DRMfldr,'RDRM.mat'),'RDRM') % uncorrected DRM
load(fullfile(DRMfldr,'DRM.mat'),'DRM') % corrected DRM

% Alter data reduction matrices if calibration wavelengths and sample
% wavelengths are different
RDRM = RDRM(:,:,calIDX);
DRM = DRM(:,:,calIDX);

% MISCELLANEOUS
NumANG = length(Angles);
NumWL = length(Wavelengths);

tol = 10^(-10); % tolerance specified as 10^(-10) at time code was passed on, Matt G.

%% CALCULATIONS
%==========================================================================
% PREALLOCATION
RMM = zeros(4,4,NumWL); 
NRMM = zeros(4,4,NumWL);
MM = zeros(4,4,NumWL); 
NMM = zeros(4,4,NumWL);
normIRR = zeros(NumWL,NumANG);

for ii = 1:NumWL
    normIRR(ii,:) = Irradiances(ii,:)/max(CALIrradiances(ii,:)); % normalize sample intensities by calibration intensities
    
    % UNCALIBRATED NORMALIZED MUELLER MATRIX 
    RMM(:,:,ii) = reshape(RDRM(:,:,ii)*normIRR(ii,:)',[4,4])'; % transposed RMM, 17Feb21 Matt Goforth/Adam Smith
    NRMM(:,:,ii) = RMM(:,:,ii)/RMM(1,1,ii); % normalize RMM by rm11
    
    % CALIBRATED NORMALIZED MUELLER MATRIX 
    MM(:,:,ii) = reshape(DRM(:,:,ii)*normIRR(ii,:)',[4,4])'; % transposed MM, 17Feb21 Matt Goforth/Adam Smith
    NMM(:,:,ii) = MM(:,:,ii)/MM(1,1,ii); % normalize MM by m11
    
    for jj = 1:4
        for kk = 1:4
            if abs(NMM(jj,kk,ii)) < 10^(-10)
                NMM(jj,kk,ii) = 0;
            elseif abs(1-NMM(jj,kk,ii)) < 10^(-10)
                NMM(jj,kk,ii) = 1;
            end
        end
    end    
end

if TOG.Comp == 1
    % DIFFERENCE OF NORMALIZED PHYSICAL MUELLER MATRIX FROM NORMALIZED IDEAL
    NMMdiff = NMM - PARA.NMMideal; 
    NMMdiff(1,1,:) = MM(1,1,:) - PARA.MMideal(1,1,:);
end

%% PLOT
%==========================================================================
if TOG.Plot == 1
    fWL = Wavelengths(1); % first element of 'VARIABLE'
    eWL = Wavelengths(end); % end element of 'VARIABLE'
    key = reshape(1:16,[4 4])'; % key for placing Mueller element in appropriate subplot
                                % position given indices 'ii' & 'jj'        

    % PLOT UNCALIBRATED NORMALIZED MUELLER MATRIX
    fig1 = figure();
    TITLE1 = 'Uncalibrated Normalized Mueller Matrix';
    PlotMueller(Wavelengths,RMM,NRMM,TITLE1,fig1,PARA,TOG.Comp,0)
   
    % PLOT NORMALIZED MUELLER MATRIX
    fig2 = figure();
    TITLE2 = 'Calibrated Normalized Mueller Matrix';
    PlotMueller(Wavelengths,MM,NMM,TITLE2,fig2,PARA,TOG.Comp,0)
    
    if TOG.Comp == 1
        % PLOT DIFFERENCE OF NMM FROM NORMALIZED IDEAL MM
        fig3 = figure();
        TITLE3 = 'Difference between Measured NMM and Ideal NMM';
        PlotMueller(Wavelengths,MM,NMMdiff,TITLE3,fig2,PARA,TOG.Comp,0)
        
        for ii = 1:4
            for jj = 1:4
                subplot(4,4,key(ii,jj))
                plot(Wavelengths,squeeze(NMMdiff(ii,jj,:)),'r','LineWidth',1)
                if ii == 1 && jj == 1
                    ylabel(sprintf('m_{%u%u} - m^{Ideal}_{%u%u}',ii,jj,ii,jj),'FontSize',12)
                else
                    ylabel(sprintf('m_{%u%u}/m_{11} - m^{Ideal}_{%u%u}/m^{Ideal}_{11}',ii,jj,ii,jj),'FontSize',12)
                end
                xlabel('\lambda (nm)','FontSize',12)
                title(sprintf('m_{%u%u}',ii,jj),'FontSize',12)
                xlim([fWL eWL])
                grid on
            end
        end
        sgtitle(sprintf('Difference between Measured NMM and Ideal NMM,  %s',PARA.SAMPLE),'FontSize',20)
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
        %==================================================================
    end
end

%% SAVE
%==========================================================================
flpthRMM = fullfile(fldr,'RMM'); 
save(flpthRMM,'RMM'); % save RMM.mat data file

flpthNRMM = fullfile(fldr,'NRMM'); 
save(flpthNRMM,'NRMM'); % save NMM.mat data file

flpthMM = fullfile(fldr,'MM'); 
save(flpthMM,'MM'); % save MM.mat data file

flpthNMM = fullfile(fldr,'NMM'); 
save(flpthNMM,'NMM'); % save NMM.mat data file

if TOG.Plot == 1
    saveas(fig1,flpthNRMM,'jpg'); % save NRMM plot as .jpg
    saveas(fig2,flpthNMM,'jpg'); % save NMM plot as .jpg
end

if TOG.Comp == 1 
    flpthNMMdiff = fullfile(fldr,'NMMdiff'); 
    save(flpthNMMdiff,'NMMdiff'); % save NMMdiff.mat data file
    if TOG.Plot == 1
        saveas(fig3,flpthNMMdiff,'jpg'); % save NMMdiff plot as .jpg
    end
end
end
