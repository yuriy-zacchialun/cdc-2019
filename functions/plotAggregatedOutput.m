function plotAggregatedOutput(T,Ts,multiy,ym,yM,ya,fullscreen,name2save)
    % % % ----------------------------------------------------------------------------------------------------
    % % Plot the system's output
    % % % ----------------------------------------------------------------------------------------------------
    % %     Inputs:
    % %         T:          time horizon for plotting
    % %         Ts:         sampling time
    % %         y:          system's output, for T+1 time steps
    % %         yNL:        system's output in noiseless setting, for T+1 time steps
    % %         ack:        sequence of acknowledgements, indicates whether control packet was received
    % %         fullscreen: equals to 1 for the fullscreen figures, otherwise 0
    % %         name2save:  string containing the desired file name and path 
    % % % ----------------------------------------------------------------------------------------------------
    %multiym = zeros(n_y,T+1,numSimulations);
    [~,~,numSims] = size(multiy);
    const = 1;
    [fSize, ~, lWidth] = getPlotSetting();
    y1M = round(max(yM(1,:)),1)+.1;
    y1m = round(min(ym(1,:)),1)-.1;
    y2M = round(max(yM(2,:)),1)+.1;
    y2m = round(min(ym(2,:)),1)-.1;
    plotName = 'Plant state on a lossy channel';
    k = 1 : T+1;          % % variable for plots involving time series
    if fullscreen == 1
        figure('Name',plotName,'units','normalized','outerposition',[0 0 1 1])
    else
        figure('Name',plotName)
    end
    ax1 = subplot(2,1,1); % noise component 1
    plot(ax1,Ts*k,multiy(1,k,1),'LineStyle','-','Color','y','LineWidth',lWidth)
    grid on
    hold all
    if numSims > 1
        for i = 2 : numSims
            plot(ax1,Ts*k,multiy(1,k,i),'LineStyle','-','Color','y','LineWidth',lWidth,'HandleVisibility','off')
        end
    end    
    plot(ax1,Ts*k,yM(1,k),'LineStyle','-','Color','b','LineWidth',lWidth)
    grid on
    % hold(ax1,'on')
    hold all
    ax1.FontSize = fSize;
    ylabel(ax1, 'Position [m]','FontSize',fSize)
    xlabel(ax1, 'Time [s]','FontSize',fSize)
    xlim(ax1,[0 const*T*Ts])
    ylim(ax1,[y1m y1M])
    title(ax1,'Cart''s position [m]','FontSize',fSize)
    plot(ax1,Ts*k,ym(1,k),'LineStyle','-','Color','g','LineWidth',lWidth)
    plot(ax1,Ts*k,ya(1,k),'LineStyle','-','Color','r','LineWidth',lWidth)
    hl1 = legend(ax1, '$\mathrm{x}_{aggreg.}$','$\mathrm{x}_{\max}$','$\mathrm{x}_{\min}$', ...
        '$\mathrm{x}_{\mathrm{avg}}$');
%     hl1 = legend(ax1, '$\mathrm{x}_{aggreg.}$','$\mathrm{x}_{\max}$','$\mathrm{x}_{\min}$', ...
%         '$\mathrm{x}_{\mathrm{avg}}$','Location','southwest');
    set(hl1, 'Interpreter', 'latex');
    hold off
    % % %
    ax2 = subplot(2,1,2); % noise component 2
    plot(ax2,Ts*k,multiy(2,k,1),'LineStyle','-','Color','y','LineWidth',lWidth)
    grid on
    hold(ax2,'on')
    if numSims > 1
        for i = 2 : numSims
            plot(ax2,Ts*k,multiy(2,k,i),'LineStyle','-','Color','y','LineWidth',lWidth,'HandleVisibility','off')
        end
    end  
    ax2.FontSize = fSize;
    ylabel(ax2, 'Angle [rad]','FontSize',fSize)
    ylim(ax2, [y2m y2M])
    xlabel(ax2,'Time [s]','FontSize',fSize)
    xlim(ax2,[0 const*T*Ts])    
    title(ax2,'Pendulum''s angle from vertical [rad]','FontSize',fSize)
    plot(ax2,Ts*k,yM(2,k),'LineStyle','-','Color','b','LineWidth',lWidth)
    plot(ax2,Ts*k,ym(2,k),'LineStyle','-','Color','g','LineWidth',lWidth)
    plot(ax2,Ts*k,ya(2,k),'LineStyle','-','Color','r','LineWidth',lWidth)
    hl2 = legend(ax2, '$\phi_{aggreg.}$','$\phi_{\max}$','$\phi_{\min}$', '$\phi_{\mathrm{avg}}$');
%     hl2 = legend(ax2, '$\phi_{aggreg.}$','$\phi_{\max}$','$\phi_{\min}$', '$\phi_{\mathrm{avg}}$', ...
%         'Location','southwest');
    set(hl2, 'Interpreter', 'latex');
    hold(ax2,'off') 
    print(name2save,'-dpng')
end