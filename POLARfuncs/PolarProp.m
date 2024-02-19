function [PolarizationProps] = PolarProp(Wavelengths,dCOMP,fldr,TOG,PARA)

% Last Edited 07 June 2021
%==========================================================================
% DESCRIPTION:
%   Extracts polarization properties from decomposed Mueller matrix
% 
% REFERENCE:
%   S. Y. Lu and R. A. Chipman, “Interpretation of Mueller matrices based on polar decomposition,” J. Opt. Soc. Am. A 13, 1106–1113 (1996)
%
% EDITED BY:
%   Matt Goforth, 16 Nov 2020
%==========================================================================

%% SETUP
%==========================================================================
% ASSIGNMENT
Dvect = dCOMP.Dvect; % diattenuater vector 
Pvect = dCOMP.Pvect; % polarizance vector
Rvect = dCOMP.Rvect; % retarder vector
depolMM = dCOMP.depolMM; % depolarizer Mueller matrix

% MISCELLANEOUS
NumWL = length(Wavelengths);

%% CALCULATIONS
%==========================================================================
% POLARIZATION PROPERTIES
LinRet = sqrt(Rvect(1,:).^2+Rvect(2,:).^2); % linear retardance
CirRet = abs(Rvect(3,:)); % circular retardance
RetOri = 0.5*atand(Rvect(2,:)./Rvect(1,:)); % retarder orientation
                
LinDiat = sqrt(Dvect(1,:).^2+Dvect(2,:).^2); % linear diattenuation
CirDiat = abs(Dvect(3,:)); % circular diattenuation
DiatOri = 0.5*atand(Dvect(2,:)./Dvect(1,:)); % diattenuater orientation

LinPolz = sqrt(Pvect(1,:).^2+Pvect(2,:).^2); % linear polarizance
CirPolz = abs(Pvect(3,:)); % circular polarizance

DepInd = zeros(1,NumWL); % depolarization index
ExtRatio = zeros(1,NumWL); % extinction ratio
for kk = 1:NumWL
    DepInd(kk) = 1-sqrt((sum(sum(depolMM(:,:,kk).^2)))-1)/sqrt(3); % depolarization index
    ExtRatio(kk) = (1+sqrt(Dvect(1,kk).^2+Dvect(2,kk).^2))/(1-sqrt(Dvect(1,kk).^2+Dvect(2,kk).^2)); % extinction ratio
end

%% PLOT RESULTS
%==========================================================================
if TOG.Plot == 1
    % SETUP
    DEG = char(176); % degree symbol
    fWL = Wavelengths(1); % first wavelength
    eWL = Wavelengths(end); % end wavelength
    PolProp = {LinRet,RetOri,CirRet,LinDiat,DiatOri,CirDiat,LinPolz,CirPolz,DepInd}; 
    labels = {sprintf('Linear Retardance (%s)',DEG),sprintf('Retarder Orientation (%s)',DEG),...
              sprintf('Circular Retardance (%s)',DEG),'Linear Diattenuation',sprintf('Diattenuater Orientation (%s)',DEG),...
              'Circular Diattenuation','Linear Polarizance','Circular Polarizance','Depolarization Index'};

    % PLOT POLARIZATION PROPERTIES
    fig1 = figure();
    for jj = 1:9
        subplot(3,3,jj)
        plot(Wavelengths,PolProp{jj},'r','LineWidth',1)
        title(labels{jj},'FontSize',12)
        xlim([fWL eWL])
        xlabel('\lambda (nm)','FontSize',12)
        grid on
    end       
    sgtitle(sprintf('Polarization Properties,  %s',PARA.SAMPLE),'FontSize',20)
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
    set(fig1,'Position', get(0,'Screensize')); % maximizes figure window; looks better when saving figure

    % PLOT EXTINCTION RATIO
    fig2 = figure();
    plot(Wavelengths,log10(ExtRatio),'r','LineWidth',1);
    xlim([fWL eWL])
    title(sprintf('Extinction Ratio,  %s',PARA.SAMPLE),'FontSize',20,'FontWeight','Normal');
    ylabel('log(Extinction Ratio)','FontSize',12);
    xlabel('\lambda (nm)','FontSize',12);
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
end

%% SAVE
%==========================================================================
% OUTPUT
% Save polarization porperties under structure
PolarizationProps.LinearRetardance = LinRet;
PolarizationProps.CircularRetardance = CirRet;
PolarizationProps.RetarderOrientation = RetOri;
PolarizationProps.LinearDiattenuation = LinDiat;
PolarizationProps.CircularDiattenuation = CirDiat;
PolarizationProps.DiattenuatorOrientation = DiatOri;
PolarizationProps.LinearPolarizance = LinPolz;
PolarizationProps.CircularPolarizance = CirPolz;
PolarizationProps.DepolarizationIndex = DepInd;
PolarizationProps.ExtinctionRatio = ExtRatio;
PolarizationProps.Wavelengths = Wavelengths;
PolarizationProps.RetardanceUnits = 'degrees';
PolarizationProps.OrientationUnits = 'degrees';
PolarizationProps.WavelengthUnits = 'nm';
save(fullfile(fldr,'PolarizationProperties'),'PolarizationProps');

if TOG.Plot == 1
    % SAVE PLOTS
    flpthPolProp = fullfile(fldr,'PolarizationProperties'); 
    saveas(fig1,flpthPolProp,'jpg'); % save PolProp plot as .jpg

    flpthExtRatio = fullfile(fldr,'ExtRatio'); 
    saveas(fig2,flpthExtRatio,'jpg'); % save ExtRatio plot as .jpg
end
end
