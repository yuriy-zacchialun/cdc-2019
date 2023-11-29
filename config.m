    % % % ----------------------------------------------------------------------------------------------------    
    % % This script defines all the parameters used by main.m
    % % % ----------------------------------------------------------------------------------------------------    
    folderMain = '/Users/yuriyzacchialun/Documents/MATLAB/cdc-2019'; % % % MacOS only, also lines 103--106
    % % % ----------------------------------------------------------------------------------------------------
    P0dBm = 0;                % % transmitter power [dBm], see [802.15.4-2006-std]. Range: from -32 to 10.
    P1dBm = 10;               % % interferer power  [dBm], see [802.15.4-2006-std]. Range: from -32 to 10.
                              % % devices are expected to operate with transmit powers between 
                              % %   -3 dBm and 10 dBm, with 0 dBm being typical: [802.15.4-2006-std], pg. 263.
                              % % % WirelessHART specification: phyTransmitPower: [BS EN 62591-2016], pg. 373:
                              % % % The values from -10 dBm to +10 dBm are supported
    % % % ----------------------------------------------------------------------------------------------------
    % % Power control parameters:
    % % % ----------------------------------------------------------------------------------------------------
    powerControl = 1;   % % 0 if the power control is disabled 
                        % % 1 if the power control is considered for both transmitters
                        % % 2 if the power control is considered only for transmitter of interest
                        % % 3 if the power control is considered only for interfearing transmitter                        
    % % % ----------------------------------------------------------------------------------------------------
    consoleMessageOff = 1;  % % % 1 suppresses the trace of execution, 0 shows it
    decimal = 6;    % % defines the number of decimal digits of round(), it is useful to present the results
    % Ts = 0.1;       % % sampling time: 100 ms is a minimum update period of a Publish data message in WHART
    Ts = 0.01;       % % sampling time: 100 ms is a minimum update period of a Publish data message in WHART
    % Ts = 0.005;       % % sampling time: 5 ms is a minimum update period of a Publish data message in WHART
    laggingAnalysis = 0;    % % If 1, the lagging effect associated with a zero-order hold discretisation is 
                            % % analyzed, 0 otherwise.
    % % % ----------------------------------------------------------------------------------------------------
    % % If we consider the inverted pendulum (on a cart) described in pendulum.m , the parameters are
    % % % ----------------------------------------------------------------------------------------------------
    M =  .5;    % % % (M)       mass of the cart                         [kg]
    m =  .2;    % % % (m)       mass of the pendulum                     [kg]
    b =  .1;    % % % (b)       coefficient of friction for cart         [N/m/sec]
    I =  .006;  % % % (I)       mass moment of inertia of the pendulum   [kg.m^2]
    g = 9.81;   % % % (g)       ravitational acceleration                [m/s^2]
    l =  .3;    % % % (l)       length to pendulum center of mass        [m]
    % % % ----------------------------------------------------------------------------------------------------
    x_0 = [0;0;pi/10;0];  % % initial state: 
                          % %   cart is motionless in the coordinate's origin and 
                          % %   pendulum lies on a cart, i.e. angle \theta is pi\2
    % % % ----------------------------------------------------------------------------------------------------
    % % The LQR problem requires definition of matrices R, Q, Sigma_w
    % % % ----------------------------------------------------------------------------------------------------
    R = 1;                          % % i.e. U in [Schenato et al, 2007]. May be defined as U = D' * D;
    Q = diag([5000,0,100,0],0);      % % i.e. W in [Schenato et al, 2007].
    Nbar = 1; % % -61.55;  % % Feedforward scaling factor 
    q = [.030;.100;.010;.150];
    Sigma_w = q * q';           % % i.e. Q in [Schenato et al, 2007].
    % % % ----------------------------------------------------------------------------------------------------
    T = 1.2e3;          % % time horizon for simulation and plotting
    % % % ---------------------------------------------------------------------------------------------------- 
    fullscreen = 1;             % % set the variable to 1 for the fullscreen figures, otherwise 0
    showPlotChannelPDF = 0;     % % set 1 to plot channel's PDF thresholds for channel/controller, otherwise 0
    % % % ----------------------------------------------------------------------------------------------------
    % % WirelessHART communication channel model parameters
    % % % ----------------------------------------------------------------------------------------------------
    % % macMaxFrameRetries = 0; % % The maximum number of retries allowed after a transmission failure, see
                            % % [802.15.4-2006-std], pg. 164, Table 86. Type: Integer, Range: from 0 to 7. 
    % % RxSensitivity = -85;    % % [dBm] See [802.15.4-2006-std], pg. 49.
    symFreq0 = 6.25e4;      % % symbol frequency [symbols/s], must be 62.5 ksymbols/s, see [802.15.4-2006-std].
    % % % fprintf('The Ts'' bandwidth is %8.3f\n\n',pow2db(6.25e4))
    W = 2e6;                % % channel bandwidth [Hz], must be 2 MHz, see [802.15.4-2006-std].
    bits4sym = 4;       % % % (16-ary ortogonal modulation such as O-QPSK encodes 4 bits with 1 symbol)
    sym4Byte = 2;       % % % (16-ary ortogonal modulation such as O-QPSK requires 2 symbols for 1 Byte)
    % % % ----------------------------------------------------------------------------------------------------
    % % Shadowing parameters [under log-normal model]:
    % % % ----------------------------------------------------------------------------------------------------
    sb0 = db2pow(2);  % % Most empirical studies support a standard deviation [dB] 
                      % % ranging from 4 to 13 dB, see [Goldsmith, 2005] at pg. 47, for OUTDOOR only. 
                      % % [yao1996] uses the standard deviation of lognormal shadowing = 2 dB.
    sb1 = sb0;        % %     
    dc0 = 8;        % % decorrelation decay distance, used to describe the correlation properties [m]
    dc1 = 20;       % % decorrelation decay distance, used to describe the correlation properties [m]
    V0 = 5.37;      % % [m/s] % % V0 = 1;
    V1 = 5.37;      % % [m/s] % % V1 = 1;
    % % % ----------------------------------------------------------------------------------------------------
    % % It is useful to remember the -174 dBm/Hz value, sometimes called excess noise ratio. It is the output 
    % % power of a 50 Ohm resistor at the average temperature of the earth (290K, 17°C) in a 1 Hz bandwidth.
    % % % ----------------------------------------------------------------------------------------------------
    k_B = 1.3806485279e-23;     % % [J/K] Boltzmann's constant
    T_ref = 290;                % % [K] temperature
    n_f = db2pow(23.8);         % % noise figure 23.8 instead of 5
    N0 = k_B * T_ref * n_f;     % % [W/Hz] noise spectral density
                                % % In case of the ideal receiver n_f = 1
    % % % ----------------------------------------------------------------------------------------------------
    te0 = 1.52e-3;  % % decorrelation time, within which the power control error is significant, see 
                    % %     [Fischione et al, 2007]
    te1 = te0;      % %
    resExp =  0.15;             % % residual power control error process' exponent
    if powerControl == 0        % % the power control is disabled 
        se0 = 0; 
        se1 = se0;
    elseif powerControl == 1    % %  the power control is considered for both transmitters
        se0 = 10^resExp;          % % standard deviation of the zero mean Gaussian process representing the 
                                % %     residual power control error process in log units
        se1 = se0;
    elseif powerControl == 2    % % the power control is considered only for transmitter of interest
        se0 = 10^resExp;
        se1 = 0;
    elseif powerControl == 3    % %  the power control is considered only for interfearing transmitter 
        se0 = 0;
        se1 = 10^resExp;    
    end
    % % % ----------------------------------------------------------------------------------------------------
    nameOutputBernoulli = 'outputBernoulli';
    nameOutputMarkov    = 'outputMarkov';
    nameOutputMarkovOracle = 'outputMarkovOracle';
    nameOutputLossless  = 'outputLossless';
    nameOutputStandard  = 'outputStandard';
    nameStatesBernoulli = 'statesBernoulli';
    nameStatesMarkov    = 'statesMarkov';
    nameStatesLossless  = 'statesLossless';
    nameStatesStandard  = 'statesStandard';
    nameMarkovChannel   = 'MarkovChannelStates';
    nameMarkovControl   = 'MarkovControlStates';
    nameCostsComparison = 'costsComparison';
    nameSnirPdf         = 'channelSnirPdf';
    % % % ----------------------------------------------------------------------------------------------------
    P0 = 1e-3 * db2pow(P0dBm);
    P1 = 1e-3 * db2pow(P1dBm);
    % % % ----------------------------------------------------------------------------------------------------        
    folderFigures =   [folderMain '/figures'];
    folderData =      [folderMain '/data'];
    folderFunctions = [folderMain '/functions'];
    folderLog =       [folderMain '/log'];    
    % % % ----------------------------------------------------------------------------------------------------
    clear 'P0dBm' 'P1dBm' 'k_B' 'T_ref' 'n_f'
    % % % ----------------------------------------------------------------------------------------------------
    cd(folderData)
    save('config.mat')
    clear 'nameOutputBernoulli' 'nameOutputMarkov' 'nameOutputLossless' 'nameOutputStandard'
    clear 'nameStatesBernoulli' 'nameStatesMarkov' 'nameStatesLossless' 'nameStatesStandard'
    clear 'nameMarkovChannel' 'nameMarkovControl' 'nameSnirPdf' 'nameCostsComparison' 'nameOutputMarkovOracle'
    % % % ----------------------------------------------------------------------------------------------------
    cd(folderLog)
    diary(['wncs-log-' datestr(now,30) '.txt'])
    cd(folderFunctions)    
    % % % ----------------------------------------------------------------------------------------------------
    % % References:
    % % % ----------------------------------------------------------------------------------------------------
    % %     [Schenato et al, 2007] Schenato, Luca, and Sinopoli, Bruno, and Franceschetti, Massimo, and 
    % %         Poolla, Kameshwar, and Sastry, S. Shankar. 
    % %             Foundations of Control and Estimation Over Lossy Networks.
    % %                 Proceedings of the IEEE, Vol. 95, No. 1, January 2007. (pg. 171)
    % % % ----------------------------------------------------------------------------------------------------
    % % For the sampling time we refer to the
    % %     [BS EN 62591-2016] BSI Standards Publication.
    % %         Industrial communication networks - Wireless communication network and communication profiles
    % %             - WirelessHART^TM. (pp. 248, 426)
    % % % ----------------------------------------------------------------------------------------------------
    % % For the characterisation of the MAC and PHY of the IEEE 802.15.4, refer to
    % %     [802.15.4-2006-std] IEEE Standard for Information technology - 
    % %         Telecommunications and information exchange between systems - 
    % %         Local and metropolitan area networks - Specific requirements. 
    % %         Part 15.4: Wireless Medium Access Control (MAC) and Physical Layer (PHY) 
    % %         Specifications for Low-Rate Wireless Personal Area Networks (WPANs).
    % % % ----------------------------------------------------------------------------------------------------
    % % For the system parameters for the inverted pendulum (on a cart), see 
    % %     [Franklin et al, 2009] G. F. Franklin, J. D. Powell, and A. Emami-Naeini. 
    % %        Feedback control of dynamic systems. Prentice Hall, 6th edition, 2009.
    % %     A useful url is
    % %         http://ctms.engin.umich.edu/CTMS/index.php?example=InvertedPendulum&section=SystemModeling
    % % % ----------------------------------------------------------------------------------------------------
    % % For the general wireless communications parameters, see
    % %     [Goldsmith, 2005] Goldsmith, Andrea: Wireless Communications. Cambridge University Press (2005)
    % % % ----------------------------------------------------------------------------------------------------
    % %     [Fischione et al, 2007] Fischione, Carlo, and Graziosi, Fabio, and Santucci, Fortunato. 
    % %         Approximation for a sum of on-off log-normal processes with wireless applications. 
    % %             IEEE Transactions on Communications, 55(9), pp. 1822–-1822 (2007)
    % % % ----------------------------------------------------------------------------------------------------
    % %     [yao1996] Yao, Shee and Geraniotis, Evaggelos. 
    % %         Optimal power control law for multimedia multirate CDMA systems.
    % %             Proceedings of Vehicular Technology Conference (VTC), vol. 1, pp. 392--396 (1996)
    % % % ----------------------------------------------------------------------------------------------------
