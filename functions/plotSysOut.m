function ...
    plotSysOut(T,Ts,multiy,ym,yM,ya,plantType,fullscreen,name2save,printOn)
% The function plots the system's output
% 
% Inputs:
%   T: the simulation time horizon;
%   Ts: the plant sampling time;
%   multiy: system outputs for each time step and each of T+1 simulations;
%   ym: minimal values of system outputs for each time step;
%   yM: maximal values of system outputs for each time step;
%   ya: average values of system outputs for each time step;
%   plantType: 'RotaryDoublePendulum' or 'LinearSinglePendulumFlexibleJoint';
%   fullscreen: logical value for the fullscreen figure plotting;
%   name2save: string containing the desired file name;
%   printOn: logical value for saving the figure in pdf.

%% 
numSims = size(multiy,3);
const = 1;
[fSize, ~, lWidth] = getPlotSetting();
y1M = round(max(yM(1,:)),1)+.1;
y1m = round(min(ym(1,:)),1)-.1;
y2M = round(max(yM(2,:)),1)+.1;
y2m = round(min(ym(2,:)),1)-.1;   
plotName = 'Plant state on a lossy channel';
k = 1 : T+1;  % variable for plots involving time series
if fullscreen == 1
    figure('Name',plotName,'units','normalized','outerposition',[0 0 1 1]);
else
    f = figure('Name',plotName,'Units','centimeters');
    set(f,'PaperSize',[20 28])
    set(f,'PaperUnits','centimeters')
    set(f,'PaperPosition',[0 0 20 28])
    set(f,'Position', [0 0 20 28])
end
if strcmp ( plantType, 'RotaryDoublePendulum' )
    y3M = round(max(yM(3,:)),1)+.1;
    y3m = round(min(ym(3,:)),1)-.1; 
    ax1 = subplot(3,1,1);               
    plot(ax1,Ts*k,multiy(1,k,1),...
        'LineStyle','-','Color','y','LineWidth',lWidth)
    grid on
    hold all
    if numSims > 1
        for i = 2 : numSims
            plot(ax1,Ts*k,multiy(1,k,i),'LineStyle','-','Color','y',...
                'LineWidth',lWidth,'HandleVisibility','off')
        end
    end
    plot(ax1,Ts*k,yM(1,k),'LineStyle','-','Color','b','LineWidth',lWidth)
    grid on
    % hold(ax1,'on')
    hold all
    ax1.FontSize = fSize;
    ylabel(ax1, 'Angle [rad]','FontSize',fSize)
    % xlabel(ax1, 'Time [s]','FontSize',fSize)
    xlim(ax1,[0 const*T*Ts])
    ylim(ax1,[y1m y1M])
    title(ax1,'Rotary arm angle [rad]','FontSize',fSize)
    plot(ax1,Ts*k,ym(1,k),'LineStyle','-','Color','g','LineWidth',lWidth)
    plot(ax1,Ts*k,ya(1,k),'LineStyle','-','Color','r','LineWidth',lWidth)
    hl1 = legend(ax1, '$\mathrm{x}_{aggreg.}$','$\mathrm{x}_{\max}$',...
        '$\mathrm{x}_{\min}$', '$\mathrm{x}_{\mathrm{avg}}$',...
        'Location','southeast');
    set(hl1, 'Interpreter', 'latex');
    hold off
    ax2 = subplot(3,1,2); % noise component 2
    plot(ax2,Ts*k,multiy(2,k,1),'LineStyle','-','Color','y',...
        'LineWidth',lWidth)
    grid on
    hold(ax2,'on')
    if numSims > 1
        for i = 2 : numSims
            plot(ax2,Ts*k,multiy(2,k,i),'LineStyle','-','Color','y',...
                'LineWidth',lWidth,'HandleVisibility','off')
        end
    end
    ax2.FontSize = fSize;
    ylabel(ax2, 'Angle [rad]','FontSize',fSize)
    ylim(ax2, [y2m y2M])
    % xlabel(ax2,'Time [s]','FontSize',fSize)
    xlim(ax2,[0 const*T*Ts])
    title(ax2,'Bottom pendulum angle from vertical [rad]','FontSize',fSize)
    plot(ax2,Ts*k,yM(2,k),'LineStyle','-','Color','b','LineWidth',lWidth)
    plot(ax2,Ts*k,ym(2,k),'LineStyle','-','Color','g','LineWidth',lWidth)
    plot(ax2,Ts*k,ya(2,k),'LineStyle','-','Color','r','LineWidth',lWidth)
    hl2 = legend(ax2, '$\phi_{aggreg.}$','$\phi_{\max}$',...
        '$\phi_{\min}$', '$\phi_{\mathrm{avg}}$','Location','northeast');
    set(hl2, 'Interpreter', 'latex');
    hold(ax2,'off')
    ax3 = subplot(3,1,3); 
    plot(ax3,Ts*k,multiy(3,k,1),'LineStyle','-','Color','y',...
        'LineWidth',lWidth)
    grid on
    hold(ax3,'on')
    if numSims > 1
        for i = 2 : numSims
            plot(ax3,Ts*k,multiy(3,k,i),'LineStyle','-','Color','y',...
                'LineWidth',lWidth,'HandleVisibility','off')
        end
    end
    ax3.FontSize = fSize;
    ylabel(ax3, 'Angle [rad]','FontSize',fSize)
    ylim(ax3, [y3m y3M])
    xlabel(ax3,'Time [s]','FontSize',fSize)
    xlim(ax3,[0 const*T*Ts])
    title(ax3,'Medium top pendulum angle from vertical [rad]',...
        'FontSize',fSize)
    plot(ax3,Ts*k,yM(3,k),'LineStyle','-','Color','b','LineWidth',lWidth)
    plot(ax3,Ts*k,ym(3,k),'LineStyle','-','Color','g','LineWidth',lWidth)
    plot(ax3,Ts*k,ya(3,k),'LineStyle','-','Color','r','LineWidth',lWidth)
    hl3 = legend(ax3, '$\phi_{aggreg.}$','$\phi_{\max}$',...
        '$\phi_{\min}$', '$\phi_{\mathrm{avg}}$');
    set(hl3, 'Interpreter', 'latex');
    hold(ax3,'off')
end
if printOn
    print(name2save,'-dpdf')
end
