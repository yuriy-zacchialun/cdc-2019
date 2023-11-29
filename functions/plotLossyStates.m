function plotLossyStates(T,Ts,x,ack,fullscreen,name2save)
    % % % ----------------------------------------------------------------------------------------------------
    % % Plot the inverted pendulum's unobservable states (velocity and angular velocity)
    % %     Inputs:
    % %         T:          time horizon for plotting
    % %         Ts:         sampling time
    % %         x:          pendulum's states, for T+1 time steps
    % %         ack:        sequence of acknowledgements, indicates whether control packet was received
    % %         fullscreen: equals to 1 for the fullscreen figures, otherwise 0
    % %         name2save:  string containing the desired file name and path 
    % % % ----------------------------------------------------------------------------------------------------
    plotName = 'Plant''s hidden states on a lossy channel';
    k = 1 : T+1;          % % variable for plots involving time series
    t = find(ack < 1);
    xLost_1 = x(2,ack < 1);
    xLost_2 = x(4,ack < 1);
    if fullscreen == 1
        figure('Name',plotName,'units','normalized','outerposition',[0 0 1 1])
    else
        figure('Name',plotName)
    end
    ax1 = subplot(2,1,1); % state component 1
    plot(ax1,Ts*k,x(2,k),'b-','LineWidth',1.5)
    grid on
    hold(ax1,'on')
    xlabel(ax1, 'Time [s]')
    ylabel(ax1, 'Velocity [m/s]')
    xlim(ax1,[0 1.2*T*Ts])
    title(ax1,'State x_2: Cart''s velocity [m/s]')
    stem(ax1,Ts*t,xLost_1,'Color','r','MarkerSize',3)
    hl1 = legend(ax1, '$\dot{x}$','\texttt{NACK}');
    set(hl1, 'Interpreter', 'latex');
    hold(ax1,'off')
    % 
    ax2 = subplot(2,1,2); % state component 2
    plot(ax2,Ts*k,x(4,k),'b-','LineWidth',1.5)
    grid on
    hold(ax2,'on')
    xlabel(ax2, 'Time [s]')
    ylabel(ax2, 'Angular velocity [rad/s]')
    xlim(ax2,[0 1.2*T*Ts])
    title(ax2,'State x_4: Pendulum''s angular velocity [rad/s]')
    stem(ax2,Ts*t,xLost_2,'Color','r','MarkerSize',3)
    hl2 = legend(ax2, '$\dot{\phi}$','\texttt{NACK}');
    set(hl2, 'Interpreter', 'latex');
    hold(ax2,'off')
    print(name2save,'-dpng')
end