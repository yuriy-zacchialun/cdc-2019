    % % % ----------------------------------------------------------------------------------------------------
    % % The main interface for the analysis of the TCP-like optimal networked state-feedback control over
    % %     WirelessHART communication channel modelled as first order Markov chain
    % % % ----------------------------------------------------------------------------------------------------
    % % The following two commands should be executed always:
    % % % ----------------------------------------------------------------------------------------------------
    run('initialisation.m')     % % prepare the environment for execution
    cd '/Users/yuriyzacchialun/Documents/MATLAB/cdc-2019'
    run('config.m')             % % define all the relevant parameters listed in the configuration file
    % % % ----------------------------------------------------------------------------------------------------
    % % % --- Select the considered scenario. ----------------------------------------------------------------
    % % % ----------------------------------------------------------------------------------------------------
    numSimulations = 1e3;
    extNumDropouts = 1; %%% It specifies wether to use the exact formula (23)=>1 or (24)=>0 of INFOCOM paper
    ell_F = 33*8;       % % % (converted to 16-ary symbols, Type 20 Command 3 for 4 variables)
    numYears = 1;
    eps_case = 2;
    if eps_case == 1        
        epsilon_b = eps;
        epsilon_p = epsilon_b;
    elseif eps_case == 2
        epsilon_p = 1/((1/Ts)*60*60*24*365.25*numYears); % % % 1 packet error in numYears of continuous operation @10Hz
        epsilon_b = epsilon_p;
    end
    numStatesChannel    = 12;      % % number of states in the Markov channel model of WirelessHART
    PartitioningBasedOnPER = 1;
    % % % ------------------------------------------------------------------------------------------------
    % % Path loss parameters:
    d0 = 10.00;       % % distance between the transmitter - receiver couple of interest [m]
    d1 =  5.00;       % % distance between the transmitter - receiver interfering couple [m]    
    % d1 =  3.50;       % % distance between the transmitter - receiver interfering couple [m]    
    % % % ----------------------------------------------------------------------------------------------------
    % % Compute the plant's system matrices sampled at proper time, or load them, if they are already computed
    % % % ----------------------------------------------------------------------------------------------------
    definePlant = 1 ;
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
    w = generateNoise(n_x,Sigma_w,T);   % % generate the process noise vector for the time horizon T+1
    cd(folderData)
    save('noise.mat','w')
    cd(folderFunctions)
    % % % -- OR ----------------------------------------------------------------------------------------------
    cd(folderData)
    % load('noise.mat')                % % load the process noise vector, if it was generated previously
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
        getAnalyticChannelSinglePersistentInterferer(...
            d0,d1,dc0,dc1,sb0,sb1,V0,V1,N0,P0,P1,symFreq0,W,te0,te1,se0,se1,Ts,ell_F,bits4sym);
    % % % [nu_Gamma,sigma_Gamma,autocov1_Gamma,autocov0der2_Gamma,autocovT1_Gamma] = ...
    % % %     getAnalyticChannelForOneInterferer(...
    % % %         d0,d1,dc,sb0,sb1,V0,V1,N0,P0,P1,symFreq0,W,te0,te1,se0,se1,Ts,ell_F);
    [expValPER,variancePER,maxNumConsecutiveDropouts,valGammaEpsp] = ...
        getAnalyticLinkQualityMetrics(nu_Gamma,sigma_Gamma,ell_F,epsilon_p,epsilon_b,...
            Ts,symFreq0,N0,W,P0,P1,se0,sb0,se1,sb1,te0,te1,V0,V1,dc0,dc1,d0,d1,extNumDropouts,bits4sym);           
    if numStatesChannel == 2 %% Gilbert Eliot channel
        GammaThresholds = valGammaEpsp;
    elseif PartitioningBasedOnPER == 0
        GammaThresholds = getEquiprobableThresholdsSNIR(nu_Gamma,sigma_Gamma,numStatesChannel);
    elseif PartitioningBasedOnPER == 1 && numStatesChannel == 3
        GammaThresholds = [getRootPER(ell_F,0.5),valGammaEpsp];
    elseif PartitioningBasedOnPER == 1 && numStatesChannel == 12
        GammaThresholds = [getRootPER(ell_F,1-epsilon_p), getRootPER(ell_F,0.9), getRootPER(ell_F,0.8), ...
            getRootPER(ell_F,0.7), getRootPER(ell_F,0.6), getRootPER(ell_F,0.5), getRootPER(ell_F,0.4), ...
            getRootPER(ell_F,0.3), getRootPER(ell_F,0.2), getRootPER(ell_F,0.1), valGammaEpsp];
    else
        error('The required setting was not implemented!')
    end    
    
    GammaSteadyStateProbability  = ...
        getSteadyStateProbabilitySNIR(nu_Gamma,sigma_Gamma,numStatesChannel,GammaThresholds);
    TPM_Analytic_symbol_rate = getMarkovChannelSymbolRateAnalytic(...
        nu_Gamma,sigma_Gamma,autocov1_Gamma,GammaThresholds,GammaSteadyStateProbability);
    TPM_Analytic_packet_rate = getMarkovChannelSymbolRateAnalytic(...
        nu_Gamma,sigma_Gamma,autocovT1_Gamma,GammaThresholds,GammaSteadyStateProbability);
    PEP_Analytic = getMarkovChannelPacketErrorProbability(numStatesChannel,nu_Gamma,...
        sigma_Gamma,ell_F,GammaThresholds,GammaSteadyStateProbability,epsilon_p);
    [eta_M_Analytic,var_M_Analytic,nStarD_AnalyticNew] = getMarkovChannelLinkQualityMetrics(...
       TPM_Analytic_packet_rate,PEP_Analytic,GammaSteadyStateProbability,epsilon_b);
    % % % ----------------------------------------------------------------------------------------------------
    ProbPacketError = PEP_Analytic;
    TPM_W = TPM_Analytic_packet_rate;
    eta_M = eta_M_Analytic;
    % % % ----------------------------------------------------------------------------------------------------
    % % Analysis of the WirelessHART communication channel
    % %     modelled as first order Markov chain with numStatesChannel states
    % % % ----------------------------------------------------------------------------------------------------
    fprintf('The probability of receiving the packet in each state of the channel''s states is\n\n\t [%5.3f', ...
        1-ProbPacketError(1))
    if length(ProbPacketError) > 2
        for i = 2 : length(ProbPacketError) - 1
            fprintf(', %5.3f', 1-ProbPacketError(i))
        end
        fprintf(', %5.3f].\n\n', 1-ProbPacketError(i+1))
    elseif length(ProbPacketError) == 2
        fprintf(', %5.3f].\n\n', 1-ProbPacketError(2))
    else
        fprintf('].\n\n')
    end
    % % % ----------------------------------------------------------------------------------------------------
    % % Test Markov channel: stabilizability and detectability analysis
    % % % ----------------------------------------------------------------------------------------------------    
    if 1 - eta_M < p_max
        warning('The average probability of receiving a packet is too low!')
    end
    % % % ----------------------------------------------------------------------------------------------------
    [A_a,B_a,C_a,P_a,doubleState,received] = generateMJLS(A,B,Q,TPM_W,ProbPacketError);
    % % % Represents the exact behaviour of the channel
    cd(folderData)
    save('sojournTime.mat','doubleState','P_a','received')
    cd(folderFunctions)
    % % % ----------------------------------------------------------------------------------------------------
    %    pause
    % OK = testMarkovChannel(A_a,B_a,C_a,P_a,consoleMessageOff);
    disp(' ')
    % p_W = p_max+eps;
    p_W = 1 - expValPER;
    [K_B, X_B] = getBernoullianController(n_x,A,B,Q,R,p_W,consoleMessageOff);
    cd(folderData)    
    save('BernoullianController.mat','K_B','X_B') % % save the outputs
    cd(folderFunctions)
    % % % ----------------------------------------------------------------------------------------------------
    % % For debug purposes only:
    [spectRadius_B, stabilizable_B] = decayConstant(A_a,B_a,P_a,K_B);
    fprintf('spectRadius_B = %11.9f\n\n', spectRadius_B)
    % % % ----------------------------------------------------------------------------------------------------
    % % Computation and analysis of LQR under Markov channel model of WirelessHART radio
    % % % ----------------------------------------------------------------------------------------------------
    % % Compute the optimal controller for Markovian packet arrivals
    % % % pause
    [K_M,J_av] = getMarkovianController(n_x,n_u,A,B,Q,R,Sigma_w,TPM_W,ProbPacketError,consoleMessageOff);
    K_b = K_M;
    for i = 1 : numStatesChannel
        K_b{i} = K_B;
    end    
    % % % ----------------------------------------------------------------------------------------------------
    cd(folderData)
    save('MarkovController.mat','K_M','J_av')
    cd(folderFunctions)
    % % % -- OR ----------------------------------------------------------------------------------------------
    % cd(folderData)
    % load('MarkovController.mat')
    % cd(folderFunctions)
    % % % ----------------------------------------------------------------------------------------------------
    disp(' ')
    [~,rho_B] =testStabilityDelModeObs(A,B,K_B,TPM_W,1-ProbPacketError);
    fprintf('SpectRadius_Lambda_delayed_B   = %5.3f\n', rho_B)
    [~,rho_M] =testStabilityDelModeObs(A,B,K_M,TPM_W,1-ProbPacketError);
    fprintf('SpectRadius_Lambda_delayed_M   = %5.3f\n\n', rho_M)
    % % % ----------------------------------------------------------------------------------------------------
    % % close all
    % % return % % % To return control to the command prompt, interrupting here the execution of the script!
    % % % ----------------------------------------------------------------------------------------------------
    
    % % % ----------------------------------------------------------------------------------------------------
    % % Markov controller behaviour for the duration of the simulation
    % % % ----------------------------------------------------------------------------------------------------
    multiack = zeros(numSimulations,T);
    multitheta = zeros(numSimulations,T+1);
    multixm = zeros(n_x,T+1,numSimulations);
    multiym = zeros(n_y,T+1,numSimulations);
    multium = zeros(n_u,T,numSimulations);
    multixb = zeros(n_x,T+1,numSimulations);
    multiyb = zeros(n_y,T+1,numSimulations);
    multiub = zeros(n_u,T,numSimulations);
    fadingLength = 0;
    maxFadingLength = 0;
    for i = 1 : numSimulations
       [ack, theta] = simulateMarkovChannel(T, numStatesChannel, TPM_W, ProbPacketError);       
       fadingLength=max(accumarray(nonzeros((cumsum(ack)+1).*~ack),1));
       if maxFadingLength < fadingLength
           maxFadingLength = fadingLength;
       end
       [x_M, y_M, u_M] = simulateGilbertStateFeedbackControl(T,A,B,C,K_M,x_0,0*w,ack,theta);
       [x_B, y_B, u_B] = simulateGilbertStateFeedbackControl(T,A,B,C,K_b,x_0,0*w,ack,theta);
       multiack(i,:) = ack;
       multitheta(i,:) = theta;
       multixm(:,:,i) = x_M;
       multiym(:,:,i) = y_M;
       multium(:,:,i) = u_M;
       multixb(:,:,i) = x_B;
       multiyb(:,:,i) = y_B;
       multiub(:,:,i) = u_B;
    end
    maxym = zeros(n_y,T+1);
    minym = zeros(n_y,T+1);
    avgym = zeros(n_y,T+1);
    maxyb = zeros(n_y,T+1);
    minyb = zeros(n_y,T+1);
    avgyb = zeros(n_y,T+1);
    for i = 1 : n_y
        for j = 1 : T+1
            maxym(i,j) = max(multiym(i,j,:));
            minym(i,j) = min(multiym(i,j,:));
            avgym(i,j) = mean(multiym(i,j,:));            
            maxyb(i,j) = max(multiyb(i,j,:));
            minyb(i,j) = min(multiyb(i,j,:));
            avgyb(i,j) = mean(multiyb(i,j,:));
        end
    end
    
    fprintf('Max fading int considered in simulations lasts ')
    fprintf('%u time steps.\n\n',maxFadingLength)
    % pause
    
    cd(folderData)
    load('config.mat','nameMarkovControl','nameStatesBernoulli')
    cd(folderFigures)
    % plotLossyOutput(T,Ts,y_M,ack,fullscreen,nameMarkovControl)
    %plotLossyOutput(T,Ts,minym,ack,fullscreen,nameMarkovControl)
    plotAggregatedOutput(T,Ts,multiym,minym,maxym,avgym,fullscreen,nameMarkovControl)
    close all
    plotAggregatedOutput(T,Ts,multiyb,minyb,maxyb,avgyb,fullscreen,nameStatesBernoulli)
    close all
    cd(folderFunctions)
    %     load('config.mat', 'nameOutputMarkov', 'nameStatesMarkov', 'nameMarkovControl')
    %     cd(folderFigures)
    %     plotLossyOutput(T,Ts,y_M,yM_NL,ack,fullscreen,nameOutputMarkov)   % % components of the plant's output vector
    % %                                                         % % i.e., position and angle, in lossy setting
    % plotLossyStates(T,Ts,x_B,ack,fullscreen,nameStatesMarkov)   % % hidden components of the plant's state,
    % %                                                         % % i.e., velocity, anglular velocity, lossy
    % plotChannelBehaviour(T,Ts,numStatesController,emme,fullscreen,nameMarkovControl)
    % % % ----------------------------------------------------------------------------------------------------
    % clear 'nameOutputMarkov' 'nameStatesMarkov' 'nameMarkovControl'
    % % % ----------------------------------------------------------------------------------------------------
    % % Performance analysis
    % % % ----------------------------------------------------------------------------------------------------
    J_B = trace(Sigma_w * X_B);
    fprintf('\t\t\tJ_B   = %-5.3f\n',J_B)
    fprintf('\t\t\tJ_M   = %-5.3f\n\n',J_av)
    % % % ----------------------------------------------------------------------------------------------------
    % % Comparison of the control gains
    % % % ----------------------------------------------------------------------------------------------------
    fprintf('    K_B   = [')
    for i = 1 : n_x
        fprintf('%5.3f, ',K_B(i))
    end
    fprintf('\b\b]\n')
    disp(' ')
    % % % ----------------------------------------------------------------------------------------------------
    for m = 1 : length(K_M)
        fprintf('    K_M(%02u) = [',m)
        for i = 1 : n_x
            fprintf('%5.3f, ',K_M{m}(i))
        end
        fprintf('\b\b]\n')
    end
    [~,~,numModes] = size(A_a);
    Q_a(:,:,numModes) = Q;
    R_a(:,:,numModes) = R;
    if numModes > 1
        for i = 1 : numModes-1
            Q_a(:,:,i) = Q;
            R_a(:,:,i) = R;
        end
    end
    % % % ----------------------------------------------------------------------------------------------------
    cd(folderMain)
    fprintf('The execution terminated on ')
    fprintf(datestr(now))
    disp(' ')
    diary OFF
