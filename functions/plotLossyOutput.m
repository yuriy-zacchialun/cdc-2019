function plotLossyOutput(T,Ts,y,ack,fullscreen,name2save)
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
    const = 1;
    [fSize, mSize, lWidth] = getPlotSetting();
    plotName = 'Plant output on a lossy channel';
    k = 1 : T+1;          % % variable for plots involving time series
    t = find(ack < 1);
    yLost_1 = y(1,ack < 1);
    yLost_2 = y(2,ack < 1);
    if fullscreen == 1
        figure('Name',plotName,'units','normalized','outerposition',[0 0 1 1])
    else
        figure('Name',plotName)
    end
    ax1 = subplot(2,1,1); % noise component 1
    plot(ax1,Ts*k,y(1,k),'LineStyle',':','Color',[0.5,0.5,0.5],'LineWidth',lWidth)
    grid on
    % hold(ax1,'on')
    hold all
    ax1.FontSize = fSize;
    ylabel(ax1, 'Position [m]','FontSize',fSize)
    xlabel(ax1, 'Time [s]','FontSize',fSize)
    xlim(ax1,[0 const*T*Ts])
    ylim(ax1,[-3.5 3.5])
    title(ax1,'Output y_1: Cart''s position [m]','FontSize',fSize)
    stem(ax1,Ts*t,yLost_1,'Color','r','MarkerSize',mSize)
    hl1 = legend(ax1, '$\mathrm{x}$ [noisy]','\texttt{NACK}','$\mathrm{x}$ [noisless]');
    set(hl1, 'Interpreter', 'latex');
    hold off
    % hold(ax1,'off')
    % % %
    ax2 = subplot(2,1,2); % noise component 2
    plot(ax2,Ts*k,y(2,k),'LineStyle',':','Color',[0.5,0.5,0.5],'LineWidth',lWidth)
    ax2.FontSize = fSize;
    grid on
    hold(ax2,'on')
    ylabel(ax2, 'Angle [rad]','FontSize',fSize)
    ylim(ax2, [-pi/4 pi/4])
    xlabel(ax2,'Time [s]','FontSize',fSize)
    xlim(ax2,[0 const*T*Ts])    
    title(ax2,'Output y_2: Pendulum''s angle from vertical [rad]','FontSize',fSize)
    stem(ax2,Ts*t,yLost_2,'Color','r','MarkerSize',mSize)
    hl2 = legend(ax2, '$\phi$ [noisy]','\texttt{NACK}', '$\phi$ [noisless]');
    set(hl2, 'Interpreter', 'latex');
    hold(ax2,'off')    
    print(name2save,'-dpng')
end