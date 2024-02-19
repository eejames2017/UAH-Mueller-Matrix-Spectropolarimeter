function [CalData] = Calibration(Angles,Wavelengths,Irradiances,fldr,TOG,PARA)

% Last Edited 07 June 2021
%==========================================================================
% DESCRIPTION:
%   Processes raw data in calibration file, calculates polarization element orientation
%   errors from Fourier decomposition of irradiance modulation through
%   least-square fit
% 
% REFERENCE:
%  P. Raman, "Spectropolarimetric characterization of light scattering materials,” 
% Ph.D. thesis, University of Alabama in Huntsville (2012): Chapter 3
%
% EDITED BY:
%   Matt Goforth, 17 Nov 2020
%==========================================================================

%% SETUP
%==========================================================================
% MISCELLANEOUS
NumWL = length(Wavelengths); % number of wavelengths scanned over by polarimeter
NumANG = length(Angles); % number of angles assumed by R1 & R2

%% MISCELLANEOUS CALCULATIONS
%==========================================================================
% IDIFF
IDIFF = 100*(Irradiances(:,1)-Irradiances(:,end))./(Irradiances(:,1)+Irradiances(:,end))/2; % difference between first and last...
                                                                                         ...measured irradiance for a given wavelength

% NORMALIZED AVERAGE POLARMETRIC IRRADIANCE MODULATION
% Measured PIM
IMOD = zeros(1,NumANG); % preallocate
for ii = 1:NumANG
    IMOD(ii) = mean(Irradiances(:,ii));
end
NIMOD = IMOD/max(IMOD);

% Ideal PIM
IMODideal61 = zeros(1,NumANG);
IMODideal361 = zeros(1,361);
for ii = 1:NumANG
    if TOG.Mode == 1
        preS = zeros(4,NumWL);
        [n,k] = RefIDXcatalog(PARA.REFIDX,Wavelengths);
        for jj = 1:NumWL
            preS(:,jj) = MMpol(0)*MMret(90,5*Angles(ii))*MMmirror(n(jj),k(jj),PARA.INCIDENTnum)*MMret(90,Angles(ii))*[1 1 0 0]';
        end
        S = mean(preS,2);
    else
        S = MMpol(0)*MMret(90,5*Angles(ii))*MMair()*MMret(90,Angles(ii))*[1 1 0 0]';
    end
    IMODideal61(ii) = S(1);
end

for ii = 1:361
    if TOG.Mode == 1
        preS = zeros(4,NumWL);
        [n,k] = RefIDXcatalog(PARA.REFIDX,Wavelengths);
        for jj = 1:NumWL
            preS(:,jj) = MMpol(0)*MMret(90,5*(ii-1))*MMmirror(n(jj),k(jj),PARA.INCIDENTnum)*MMret(90,(ii-1))*[1 1 0 0]';
        end
        S = mean(preS,2);
    else
        S = MMpol(0)*MMret(90,5*(ii-1))*MMair()*MMret(90,(ii-1))*[1 1 0 0]';
    end
    IMODideal361(ii) = S(1);
end
NIMODideal61 = IMODideal61/max(IMODideal61); 
NIMODideal361 = IMODideal361/max(IMODideal361);

%% LEAST-SQAURE FIT
%==========================================================================
% Calculate polarimeter's systemic errors through least-squares fit.

if TOG.Mode == 1
    % For 'ModelFunctionREFL'
    Cal0 = [0 0 0 0 0 0 0 0]; % initialization coefficients for 'lsqcurvefit'                  
    ub =   []; % upper bounds for solution given by 'lsqcurvefit'
    lb =   []; % lower bounds for solution given by 'lsqcurvefit'
    % [R1 retardance error, R2 retardance error, R1 orientation error, R2 orientation error, P2 orientation error,...
    ...MIR real refractive index error, MIR imaginary refractive index error, MIR orientation error]

    CalData = zeros(NumWL,8);
    OPTIONS = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt'); % added by Matt G., 22 Oct 2020
    [n,k] = RefIDXcatalog(PARA.REFIDX,Wavelengths); 
    for ii = 1:NumWL
        normIRR = Irradiances(ii,:)/max(Irradiances(ii,:)); % vector of irradiances for each R1 angle; ...
                                                              ...each row normalized by its largest valued irradiance 
        x0 = {Angles,n(ii),k(ii),PARA.INCIDENTnum};
        CalData(ii,1:end) = lsqcurvefit(@ModelFunctionREFL, Cal0, x0, normIRR, lb, ub, OPTIONS);
    end
else
    % For 'ModelFunctionTRANS'
    Cal0 = [0 0 0 0 0]; % initialization coefficients for 'lsqcurvefit'                  
    ub =   []; % upper bounds for solution given by 'lsqcurvefit'
    lb =   []; % lower bounds for solution given by 'lsqcurvefit'
    % [R1 retardance error, R2 retardance error, R1 orientation error, R2 orientation error, P2 orientation error]

    CalData = zeros(NumWL,5);
    OPTIONS = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt'); % added by Matt G., 22 Oct 2020
    for ii = 1:NumWL
        normIRR = Irradiances(ii,:)/max(Irradiances(ii,:)); % vector of irradiances for each R1 angle; ...
                                                              ...each row normalized by its largest valued irradiance 
        CalData(ii,1:end) = lsqcurvefit(@ModelFunctionTRANS, Cal0, Angles, normIRR, lb, ub, OPTIONS);
    end
end

%% MAKE PLOTS / ANNOTATE
%==========================================================================
if TOG.Plot == 1
    % SETUP
    DEG = char(176); % degree symbol
    labels = {sprintf('Retardance Error of R1 (%s)',DEG),sprintf('Retardance Error of R2 (%s)',DEG),...
              sprintf('Orientation Error of R1 (%s)',DEG),sprintf('Orientation Error of R2 (%s)',DEG),...
              sprintf('Orientation Error of P2 (%s)',DEG),'Real Ref. Idx. n_{Mirror} Error',...
              'Imaginary Ref. Idx. k_{Mirror} Error',sprintf('Orientation Error of MIR (%s)',DEG)};
    fWL = Wavelengths(1); % first wavelength
    eWL = Wavelengths(end); % end wavelength
    [~,NumCal] = size(CalData);

    % PLOT IMOD          
    fig1 = figure();
    plot(Angles,NIMOD,'-rx','LineWidth',1,'MarkerEdgeColor','k','MarkerSize',7,'MarkerFaceColor','k')
    hold on
    plot(Angles,NIMODideal61,'-b','LineWidth',0.5)
    plot(0:1:360,NIMODideal361,'-g','LineWidth',0.5)
    hold off
    xlim([Angles(1) Angles(end)])
    ylim([0 1])
    xticks(Angles(1):30:Angles(end))
    yticks(0:0.2:1)
    xlabel(sprintf('Retarder Orientation \\theta (R_{1}\\theta, R_{2}5\\theta) (%c)',DEG),'FontSize',12)
    ylabel('Normalized Irradiance','FontSize',12)
    title(sprintf('Normalized \\lambda-Averaged Polarimetric Irradiance Modulation,  %s',PARA.SAMPLE),'FontSize',20,'FontWeight','Normal')
    grid on
    grid minor
    ax1 = gca;

    hold on
    pa = plot(nan,nan,'ko',nan,nan,'ko',nan,nan,'ko',nan,nan,'ko',nan,nan,'ko',nan,nan,'ko'); % create 2nd legend for listing lab setup parameters
    hold off
    
    ax1POSITION = ax1.Position; % position of first axes
    ax2 = axes('position',ax1POSITION,'visible','off'); % create 2nd invisible axes for 2nd legend
    LEG1 = legend(ax1,'Measured','Ideal (61)','Ideal (361)');
    LEG1.FontSize = 12;  
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
    
    % PLOT IDIFF
    fig2 = figure();
    plot(Wavelengths,IDIFF,'-rx','LineWidth',1,'MarkerEdgeColor','k','MarkerSize',7,'MarkerFaceColor','k')
    title(sprintf('Percent Difference between First and Last Irradiance,  %s',PARA.SAMPLE),'FontSize',20,'FontWeight','Normal')
    xlim([fWL eWL])
    xlabel('\lambda (nm)','FontSize',12)
    ylabel('Percent Difference','FontSize',12)
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
    
    % PLOT ERRORS
    fig3 = figure();
    for ii = 1:NumCal
        if TOG.Mode == 1
            subplot(4,2,ii)
        else
            subplot(3,2,ii)
        end
        plot(Wavelengths,CalData(:,ii),'r','LineWidth',1)
        title(labels{ii},'FontSize',12)
        xlim([fWL eWL])
        xlabel('\lambda (nm)','FontSize',12)
        grid on
    end       
    sgtitle(sprintf('Calibration Errors,  %s',PARA.SAMPLE),'FontSize',20)
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
end

%% SAVE 
%==========================================================================
flpthPIM = fullfile(fldr,'PIM'); 
save(flpthPIM,'NIMOD'); % save NIMOD.mat data file

flpthIDIFF = fullfile(fldr,'IDIFF'); 
save(flpthIDIFF,'IDIFF'); % save NIDIFF.mat data file

flpthCalERR = fullfile(fldr,'CalData');
save(flpthCalERR,'CalData'); % save CalData.mat data file

if TOG.Plot == 1
    saveas(fig1,flpthPIM,'jpg'); % save PIM plot as .jpg
    saveas(fig2,flpthIDIFF,'jpg'); % save IDIFF plot as .jpg
    saveas(fig3,flpthCalERR,'jpg'); % save CalData plot as .jpg
end

end
