    % % % ----------------------------------------------------------------------------------------------------
    % % The main interface for the analysis of the TCP-like optimal networked state-feedback control over

    % %     WirelessHART communication channel modelled as first order Markov chain
    % % % ----------------------------------------------------------------------------------------------------
    % % The following two commands should be executed always:
    % % % ----------------------------------------------------------------------------------------------------
    run('initialisation.m')     % % prepare the environment for execution
    cd '/Users/yuriyzacchialun/Documents/MATLAB/WNCS-CDC19/WNCS'
    run('config.m')             % % define all the relevant parameters listed in the configuration file
    % % % ----------------------------------------------------------------------------------------------------
    % % % --- Select the considered scenario. ----------------------------------------------------------------
    % % % ----------------------------------------------------------------------------------------------------
    numSimulations = 1e3;
    extNumDropouts = 0;
    ell_F = 208;
    numYears = 1;
    eps_case = 2;
    if eps_case == 1        
        epsilon_b = eps;
        epsilon_p = epsilon_b;
    elseif eps_case == 2
        epsilon_p = 1/((1/Ts)*60*60*24*365.25*numYears); % % % 1 packet error in numYears of continuous operation @10Hz
        epsilon_b = epsilon_p;
    end
    numStatesChannel    = 2;      % % number of states in the Markov channel model of WirelessHART
    % % % ------------------------------------------------------------------------------------------------
    % % Path loss parameters:
    d0 = 10.00;       % % distance between the transmitter - receiver couple of interest [m]
    d1 =  3.50;       % % distance between the transmitter - receiver interfering couple [m]    
    % % % ----------------------------------------------------------------------------------------------------
    % % Compute the plant's system matrices sampled at proper time, or load them, if they are already computed
    % % % ----------------------------------------------------------------------------------------------------
    definePlant = 0 ;
    if definePlant == 1
        [A,B,C,p_max,n_x,n_u,n_y] = getPendulum(M,m,b,I,g,l,Nbar,Ts,decimal,laggingAnalysis);
        lengthBernoulli = analyzePlant(A,B,C,p_max,n_x,n_u,n_y,epsilon_b,folderData);
    else
        cd(folderData)
        load('pendulum.mat')        % % % load system matrices, if they were computed previously
        fprintf('Max fading int under which the system remains controllable with probability 1 ')
        fprintf('lasts %u time steps.\n\n',lengthBernoulli)
    end
    clear 'M' 'm' 'b' 'I' 'g' 'l'
    cd(folderFunctions)
    % % % ----------------------------------------------------------------------------------------------------
    % w = generateNoise(n_x,Sigma_w,T);   % % generate the process noise vector for the time horizon T+1
    % cd(folderData)
    % save('noise.mat','w')
    % cd(folderFunctions)
    % % % -- OR ----------------------------------------------------------------------------------------------
    cd(folderData)
    load('noise.mat')                % % load the process noise vector, if it was generated previously
    cd(folderFunctions)
    % % % ----------------------------------------------------------------------------------------------------
    % cd(folderFigures) 
    % plotNoise(T,Ts,w,fullscreen)        % % all four components of the process noise vector
    % cd(folderFunctions)
    % % % ----------------------------------------------------------------------------------------------------
    showConfigurationSummary(d0,d1,powerControl)
    % % % ----------------------------------------------------------------------------------------------------
    % % Analysis of the analytic model of the WirelessHART communication channel:
    [nu_Gamma,sigma_Gamma,autocov1_Gamma,autocov0der2_Gamma,autocovT1_Gamma] = ...
        getAnalyticChannelForOneInterferer(...
            d0,d1,dc,sb0,sb1,V0,V1,N0,P0,P1,symFreq0,W,te0,te1,se0,se1,Ts,ell_F);    
    
    % % % ----------------------------------------------------------------------------------------------------
    % % Plot presentation
    % % % ----------------------------------------------------------------------------------------------------
    
    sinr_db = -40 : 0.1 : 40;    
    sinr_pdf_db = normpdf(sinr_db,nu_Gamma,sigma_Gamma);
    sinr_cdf_db = normcdf(sinr_db,nu_Gamma,sigma_Gamma);
    R_packet = getPER(sinr_db,ell_F);
    
    
    indLatTh = find(sinr_cdf_db >= 1 - epsilon_p,1);
    
    plotName = 'Probabilistic characterization of the SINR [dB]';
    if fullscreen == 1
        figure('Name',plotName,'units','normalized','outerposition',[0 0 1 1])
    else
        figure('Name',plotName)
    end
    [fSize, mSize, lWidth] = getPlotSetting();
    % epsilon_p
    pHeader = plot(sinr_db,sinr_pdf_db,'LineStyle','-','Color','c','LineWidth',lWidth);
    hold all
    plot(sinr_db,sinr_cdf_db,'LineStyle','-','Color','b','LineWidth',lWidth);
    plot(sinr_db,R_packet,'LineStyle','-','Color','m','LineWidth',lWidth);
    grid on
    eps_p = num2str(epsilon_p,'%5.3e\n');
    hl2 = xline(sinr_db(indLatTh),'LineStyle','--','Color','k','LineWidth',lWidth,...
        'FontSize',fSize,'Label',{'PEP = $\varepsilon$',['$\varepsilon=$~' eps_p]},'LabelVerticalAlignment','middle',...
        'LabelOrientation','horizontal');
    set(hl2, 'Interpreter', 'latex');
    ax1 = gca;
    ax1.FontSize = fSize;
    ax1.XLabel.String = 'SINR [dB]';
    ax1.YLabel.String = 'Probability';
    title(ax1,'Probabilistic characterization a wireless link','FontSize',fSize)
    hl1 = legend(ax1, '{\bf Probability distribution function}',...
        '{\bf Cumulative distribution function}',...
        '{\bf Packet error probability (PEP)}',...
        'Location', 'northwest');
    set(hl1, 'Interpreter', 'latex');    
    cd(folderFigures)
    print('probablitiesWirelessLink','-dpng')    
    cd(folderMain)
    fprintf('The execution terminated on ')
    fprintf(datestr(now))
    disp(' ')
    diary OFF
    
    function y = getPER(x,ellF) % ell_F,epsilon_p
    y = 1 - (4*exp(-10*10.^(x/10)) + (182*exp(-15*10.^(x/10)))/3 - (728*exp(-16*10.^(x/10)))/5 + ...
        (4004*exp(-18*10.^(x/10)))/15 + 429*exp(-(35*10.^(x/10))/2) - (56*exp(-(40*10.^(x/10))/3))/3 + ...
        (4004*exp(-(50*10.^(x/10))/3))/15 + (182*exp(-(55*10.^(x/10))/3))/3 - ...
        (8*exp(-(56*10.^(x/10))/3))/15 + exp(-(75*10.^(x/10))/4)/30 - (1144*exp(-(120*10.^(x/10))/7))/3 + ...
        4*exp(-(130*10.^(x/10))/7) - (1144*exp(-(160*10.^(x/10))/9))/3 - ...
        (728*exp(-(200*10.^(x/10))/11))/5 - (56*exp(-(240*10.^(x/10))/13))/3 - 1).^ellF;
    end