function [nuGamma,sigmaGamma,autocovGamma1,autocov0der2,autocovGammaT1,M1] = ...
    getAnalyticChannelSinglePersistentInterferer(...
            d0,d1,dc0,dc1,sb0,sb1,V0,V1,N0,P0,P1,symFreq0,W,te0,te1,se0,se1,Ts,ellF,bits4sym)
    % % % ----------------------------------------------------------------------------------------------------
    % % Produces the first order Markov chain representation of the wireless HART channel model with
    % % Input: 
    % %         d0: distance between the transmitter - receiver couple of interest [m]
    % %         d1: distance between the transmitter - receiver interfering couple [m]
    % % % % Shadowing parameters:
    % %         dc: decorrelation decay distance, used to describe the correlation properties [m]
    % %        sb0: standard deviation of lognormal shadowing [dB], for the Tx - Rx couple of interest
    % %        sb1: standard deviation of lognormal shadowing [dB], for the Tx - Rx interfering couple
    % %         V0: device speed [m/s], for the Tx - Rx couple of interest
    % %         V1: device speed [m/s], for the Tx - Rx interfering couple
    % % % % Standard parameters:
    % %         N0: noise spectral density [W/Hz]
    % %         P0: transmitter power [W], for the transmitter of interest
    % %         P1: transmitter power [W], for the interferer
    % %   symFreq0: symbol frequency [symbols/s]
    % %          W: channel bandwidth [Hz]
    % % % % Power control (PC) parameters:
    % %        te0: decorrelation time [s], for which the PC error is significant, for the couple of interest
    % %        te1: decorrelation time [s], for which the PC error is significant, for the interfering couple
    % %        se0: std dev of the zero mean Gaussian process representing the residual PC error process for the couple of interest [dB]
    % %        se1: std dev of the zero mean Gaussian process representing the residual PC error process for the interfering couple [dB]
    % % % % Packet-level analysis requires the study of autocorrelation between two successive packets
    % %         Ts: sampling time, i.e. interval between two packet transmissions
    % %       ellF: length of frames encapsulating data packets of interest
    % % Outputs:
    % %           nu: expected value of SINR
    % %        sigma: standard deviation of SINR
    % %     autocov1: autocovariance for tau = 1 of the Gaussian process reperesnting SNIR
    % %     autocov0: autocovariance for tau = 0 of the Gaussian process reperesnting SNIR
    % % autocov0der2: second derivative of autocovariance for tau = 0 of the Gaussian process potraying SNIR
    % % % ----------------------------------------------------------------------------------------------------
    fprintf('The characteristics of the considered analytical WirelessHART channel are \n\n')
    % % % ----------------------------------------------------------------------------------------------------
    T0 = 1/symFreq0; % symbol time
    % % Path loss parameters are obtained from [802.15.4-2006-std], pp. 273--274:
    if d0 > 8
        alfa0 = db2pow(-58.5 - 33*log10(d0/8));
    elseif d0 >= 0.5
        alfa0 = db2pow(-40.2 - 20*log10(d0));
    else
        error('Path loss coefficient is not computable due to near-field and implementation effects')
    end
    if d1 > 8
        alfa1 = db2pow(-58.5 - 33*log10(d1/8));
    elseif d1 >= 0.5
        alfa1 = db2pow(-40.2 - 20*log10(d1));
    else
        error('Path loss coefficient is not computable due to near-field and implementation effects')
    end
    % % fprintf('      alpha0 = %5.2f [dB] is the path loss at d0 = %5.2f [m]\n\n', pow2db(alfa0), d0)
    % % % ----------------------------------------------------------------------------------------------------    
    D = N0/(4*P0*alfa0^2);
    G0 = W*T0;
    B = 8*P1*alfa1^2/(3*P0*alfa0^2*G0);    
    sy0 = sqrt(se0^2+sb0^2);
    sy1 = sqrt(se1^2+sb1^2);
    % % % ----------------------------------------------------------------------------------------------------
    t = 0;
    % % % ----------------------------------------------------------------------------------------------------
    cse0 = se0^2*exp(-1/2*(t/te0)^2);   % % % Warning: t = 0
    cse1 = se1^2*exp(-1/2*(t/te1)^2);   % % % Warning: t = 0
    cb0 = sb0^2*exp(-1/2*(V0*t/dc0)^2);  % % % Warning: t = 0
    cb1 = sb1^2*exp(-1/2*(V1*t/dc1)^2);  % % % Warning: t = 0
    cy0 = cse0+cb0;
    cy1=cse1+cb1;
    % % % --- M1 and M2 --------------------------------------------------------------------------------------
    M1 = D*exp(1/2*sy0^2)+B.*exp(1/2*(sy1^2+sy0^2));
    M2_0 = exp(2*sy0^2)*(D^2+2*B*D*exp(1/2*sy1^2) + B^2*exp(2*sy1^2));
    % % % --- Components of M1 and M2 ------------------------------------------------------------------------
    X = exp(sy0^2+cy0)*D^2;                          % % % Warning: t = 0
    Y = 2*B*D*exp(1/2*sy1^2+sy0^2+cy0);              % % % Warning: t = 0
    Z = B^2*exp(sy1^2+cy1+sy0^2+cy0);                % % % Warning: t = 0
    M2 = X+Y+Z;                                      % % % Warning: t = 0
    % % % --- First derivatives of the components ------------------------------------------------------------
    cse0der1 = -t/(te0^2)*cse0;                      % % % Warning: t = 0
    cse1der1 = -t/(te1^2)*cse1;                      % % % Warning: t = 0
    cb0der1 = -(V0/dc0)^2*t*cb0;                      % % % Warning: t = 0
    cb1der1 = -(V1/dc1)^2*t*cb1;                      % % % Warning: t = 0
    cy0der1 = cse0der1+cb0der1;                      % % % Warning: t = 0
    cy1der1 = cse1der1+cb1der1;                      % % % Warning: t = 0
    % % % --- Second derivatives of the components -----------------------------------------------------------
    cse0der2 = -t/(te0^2)*cse0der1 -1/te0^2*cse0;   % % % Warning: t = 0
    cse1der2 = -t/(te1^2)*cse1der1 -1/te1^2*cse1;   % % % Warning: t = 0
    cb0der2 = -(V0/dc0)^2*t*cb0der1-V0^2/dc0^2*cb0;   % % % Warning: t = 0
    cb1der2 = -(V1/dc1)^2*t*cb1der1-V1^2/dc1^2*cb1;   % % % Warning: t = 0
    cy0der2 = cse0der2+cb0der2;                     % % % Warning: t = 0
    cy1der2 = cse1der2+cb1der2;                     % % % Warning: t = 0
    % % % --- First derivative of M2 -------------------------------------------------------------------------
    X1 = X*cy0der1;                                 % % % Warning: t = 0
    Y1 = Y*cy0der1;                                 % % % Warning: t = 0
    Z1 = Z*(cy0der1+cy1der1);                       % % % Warning: t = 0
    der1M2 = X1+Y1+Z1;                              % % % Warning: t = 0
    % % % --- Second derivative of M2 ------------------------------------------------------------------------
    X2 = X1*cy0der1+X*cy0der2;                      % % % Warning: t = 0
    Y2 = Y1*cy0der1+Y*cy0der2;                      % % % Warning: t = 0
    Z2 = Z1*(cy0der1+cy1der1)+Z*(cy0der2+cy1der2);  % % % Warning: t = 0
    der2M2 =X2+Y2+Z2;                               % % % Warning: t = 0
    % % % --- Moment matching --------------------------------------------------------------------------------
    meanZm = 2*log(M1)-1/2*log(M2_0);
    stdDevZm2 = log(M2_0)-2*log(M1);
    % % % autocovZm = log(M2)-2*meanZm-stdDevZm2;         % % % Warning: t = 0
    nuGamma  = -5/log(10)*meanZm;
    fprintf('          nu = %9.3f [dB], i.e. expected value of SINR\n', nuGamma)
    % % % sigma = (5/log(10))*sqrt(autocovZm);
    sigmaGamma = (5/log(10))*sqrt(stdDevZm2);
    fprintf('       sigma = %9.3f [dB], i.e. standard deviation of SINR\n', sigmaGamma)
    fprintf('     sigma^2 = %9.3f [dB], i.e. variance of SINR\n', (5/log(10))^2*stdDevZm2)
    % % fprintf('   sigma_old = %8.3f [dB], i.e. standard deviation of SINR\n\n', (5/log(10))*sqrt(autocovZm))
    % %
    % % % --- Autocovariance of Gamma ------------------------------------------------------------------------
    % % % --------- for t = 1 -------------------------------------------------------------------------------- 
    t = 1;
    % % % ----------------------------------------------------------------------------------------------------
    sy0t = se0^2+sb0^2;
    sy1t = se1^2+sb1^2; 
    cy0t = getcse(t,se0,te0)+getcb(t,sb0,V0,dc0);
    cy1t=getcse(t,se1,te1)+getcb(t,sb1,V1,dc1);   
    Xt = exp(sy0t+cy0t)*D^2;
    Yt = 2*B*D*exp(1/2*sy1t+sy0t+cy0t);
    Zt = B^2*exp(sy1t+cy1t+sy0t+cy0t);
    M2t = Xt+Yt+Zt;
    % % % ----------------------------------------------------------------------------------------------------
    M2_1 = M2t;
    % % % --- Autocovariance of Gamma ------------------------------------------------------------------------
    % % % --------- for t = T --------------------------------------------------------------------------------
    T = bits4sym*Ts*symFreq0/ellF;
    t = T;
    % % % ----------------------------------------------------------------------------------------------------    
    sy0t = se0^2+sb0^2;
    sy1t = se1^2+sb1^2; 
    cy0t = getcse(t,se0,te0)+getcb(t,sb0,V0,dc0);
    cy1t=getcse(t,se1,te1)+getcb(t,sb1,V1,dc1);   
    Xt = exp(sy0t+cy0t)*D^2;
    Yt = 2*B*D*exp(1/2*sy1t+sy0t+cy0t);
    Zt = B^2*exp(sy1t+cy1t+sy0t+cy0t);
    M2t = Xt+Yt+Zt;
    % % % ----------------------------------------------------------------------------------------------------
    M2_t = M2t;
    % % % ----------------------------------------------------------------------------------------------------        
    autocovGamma1 = (5/log(10))^2*log(M2_1/(M1^2));
    fprintf('  autocov(1) = %9.3f [dB], i.e. autocovariance of SINR with tau = 1\n', autocovGamma1)
    autocovGammaT1 = (5/log(10))^2*log(M2_t/(M1^2));
    fprintf('  autocov(T) = %9.3f [dB], i.e. autocovariance of SINR with tau = T\n', autocovGammaT1)
    % % % ----------------------------------------------------------------------------------------------------    
    % % % --- Autocovariance of Gamma ------------------------------------------------------------------------
    % % % --------- for t = 2T --------------------------------------------------------------------------------
    t = 2*T;
    % % % ----------------------------------------------------------------------------------------------------    
    sy0t = se0^2+sb0^2;
    sy1t = se1^2+sb1^2; 
    cy0t = getcse(t,se0,te0)+getcb(t,sb0,V0,dc0);
    cy1t=getcse(t,se1,te1)+getcb(t,sb1,V1,dc1);   
    Xt = exp(sy0t+cy0t)*D^2;
    Yt = 2*B*D*exp(1/2*sy1t+sy0t+cy0t);
    Zt = B^2*exp(sy1t+cy1t+sy0t+cy0t);
    M2t = Xt+Yt+Zt;
    % % % ----------------------------------------------------------------------------------------------------
    M2_T = M2t;
    % % % ----------------------------------------------------------------------------------------------------
    autocovGammaT2 = (5/log(10))^2*log(M2_T/(M1^2));
    fprintf(' autocov(2T) = %9.3f [dB], i.e. autocovariance of SINR with tau = 2T\n', autocovGammaT2)
    % % % --- Level Crossing Rate (LCR) parameters evaluation ------------------------------------------------  
    % % % autocov0     = (5/log(10))^2*autocovZm;   % % % By definition it equals sigmaGamma^2                       
    autocov0der2 = -(5/log(10))^2*((der2M2*M2-der1M2^2)/M2^2);                            
    % % % fprintf('The value of the lambda0 is %7.4f. \n',autocov0);
    fprintf('autocov0der2 = %9.3f [dB], i.e. second derivative of autocovariance of SINR in tau = 0\n\n',...
        autocov0der2);
    fprintf('NumJumpExact = %9.3f \n\n',T);
    % % % ----------------------------------------------------------------------------------------------------
    % % % --- Auxiliary functions ----------------------------------------------------------------------------
    function cbt = getcb(t,sb,V,dc)
        cbt = sb^2*exp(-1/2*(V*t/dc)^2);
    end
    function cset = getcse(t,se,te)
        cset = se^2*exp(-1/2*(t/te)^2);        
    end
    % % % ----------------------------------------------------------------------------------------------------
end
    % % % ----------------------------------------------------------------------------------------------------
    % % For the characterisation of the MAC and PHY of the IEEE 802.15.4, refer to
    % %     [802.15.4-2006-std] IEEE Standard for Information technology - 
    % %         Telecommunications and information exchange between systems - 
    % %         Local and metropolitan area networks - Specific requirements. 
    % %         Part 15.4: Wireless Medium Access Control (MAC) and Physical Layer (PHY) 
    % %         Specifications for Low-Rate Wireless Personal Area Networks (WPANs).
    % % % ----------------------------------------------------------------------------------------------------