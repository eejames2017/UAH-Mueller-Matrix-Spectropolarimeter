function PlotMueller(VARIABLE,MM,NMM,TITLE,figHANDLE,PARA,TOGComp,TOGDiff)

% Last Edited 25 Mar 2021
%==========================================================================
% DESCRIPTION:
%   Function to plot uncalibrated, calibrated, and physical Mueller matrix
% 
% REFERENCE:
%   ...
%
% CREATED BY:
%   Matt Goforth, 11 Mar 2021
%==========================================================================

%% SETUP
%==========================================================================
fVAR = VARIABLE(1); % first element of 'VARIABLE'
eVAR = VARIABLE(end); % end element of 'VARIABLE'
key = reshape(1:16,[4 4])'; % key for placing Mueller element in appropriate subplot
                            % position given indices 'ii' & 'jj'
                            
%% PLOT
%==========================================================================
for ii = 1:4
    for jj = 1:4
        if ii == 1 && jj == 1
            subplot(4,4,key(ii,jj))
            plot(VARIABLE,squeeze(MM(ii,jj,:)),'r','LineWidth',1)
            if TOGComp == 1
                hold on
                plot(VARIABLE,squeeze(PARA.MMideal(ii,jj,:)),'--g','LineWidth',1)
                hold off              
            end
            if TOGDiff == 1
                ylabel(sprintf('m_{%u%u} - HI',ii,jj),'FontSize',12)
            else
                ylabel(sprintf('m_{%u%u}',ii,jj),'FontSize',12)
            end
        else
            subplot(4,4,key(ii,jj))
            plot(VARIABLE,squeeze(NMM(ii,jj,:)),'r','LineWidth',1);
            if TOGComp == 1
                hold on
                plot(VARIABLE,squeeze(PARA.NMMideal(ii,jj,:)),'--g','LineWidth',1)
                hold off
            end
            ylim([-1 1]);
            if TOGDiff == 1
                ylabel(sprintf('m_{%u%u}/N - HI/N',ii,jj),'FontSize',12)
            else
                ylabel(sprintf('m_{%u%u}/m_{11}',ii,jj),'FontSize',12)
            end
        end
        xlim([fVAR eVAR])
        xlabel('\lambda (nm)','FontSize',12)
        title(sprintf('m_{%u%u}',ii,jj),'FontSize',12)
        grid on
    end
end       

sgtitle(sprintf('%s,  %s',TITLE,PARA.SAMPLE),'FontSize',20)

subplot(4,4,1)
ax1 = gca;
hold on
pa = plot(nan,nan,'ko',nan,nan,'ko',nan,nan,'ko',nan,nan,'ko',nan,nan,'ko',nan,nan,'ko'); % create 2nd legend for listing lab setup parameters
hold off

ax1POSITION = ax1.Position; % position of first axes
ax2 = axes('position',ax1POSITION,'visible','off'); % create 2nd invisible axes for 2nd legend
if TOGComp == 1
    LEG1 = legend(ax1,'Measured','Ideal');
    LEG1.FontSize = 11;
end
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

set(figHANDLE,'Position', get(0,'Screensize')); % maximizes figure window; looks better when saving figure 

end

