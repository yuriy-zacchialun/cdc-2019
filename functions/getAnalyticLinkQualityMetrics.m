function [expValPER,variancePER,maxNumConsecutiveDropouts,valGammaEpsp] = ...
	getAnalyticLinkQualityMetrics(mu,sigma,ellF,epsp,epsb,...
        Ts,symFreq0,N0,W,P0,P1,se0,sb0,se1,sb1,te0,te1,V0,V1,dc0,dc1,d0,d1,...
        extNumDropouts,b4sym)
   
    dmax = 2500; % % % Stoppage criterion for getmvncdf
    
    [expValPER, variancePER] = getMomentsPER(mu,sigma,-Inf,Inf,ellF);
    valGammaEpsp = getRootPER(ellF,epsp);
    
    fprintf('valGammaEpsp = %9.3f [dB], i.e. root of SINR: PER == eps_p\n\n',valGammaEpsp);
    maxNumConsecutiveDropouts = ceil(log(epsb)/log(0.5*(1+erf((valGammaEpsp-mu)/(sigma*sqrt(2))))));
    fprintf('ell_B^starLB = %5u.          i.e. max num of consecutive dropouts; it depends on eps_b\n\n',...
            maxNumConsecutiveDropouts)
    if extNumDropouts == 1
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
        % % % ----------------------------------------------------------------------------------------------------
        D = N0/(2*P0*alfa0^2);
        G0 = W*T0;
        B = 16*P1*alfa1^2/(3*P0*alfa0^2*G0);
        sy0 = sqrt(se0^2+sb0^2);
        sy1 = sqrt(se1^2+sb1^2);
        % % % ----------------------------------------------------------------------------------------------------
        M1 = D*exp(1/2*sy0^2)+B*exp(1/2*(sy1^2+sy0^2));
        % % % ----------------------------------------------------------------------------------------------------
        % % % --- Autocovariance of Gamma ------------------------------------------------------------------------
        % % % --------- for t = k*ell_E+1 ------------------------------------------------------------------------
        % % % ----------------------------------------------------------------------------------------------------
        T1 = b4sym*Ts*symFreq0/ellF;
        M2t = getM2t(T1);
        autocovGammaT1(1) = (5/log(10))^2*log(M2t/(M1^2));
        fprintf('autocov(%3u) = %18.15f [dB], i.e. autocovariance of SINR with tau = %3u\n',1,autocovGammaT1(1),1)
        autocovGammaT1(dmax)=0;
        i = 2;
        while autocovGammaT1(i-1) > 10*eps && i < dmax
            T = i * T1;
            M2t = getM2t(T);
            autocovGammaT1(i) = (5/log(10))^2*log(M2t/(M1^2));
            fprintf('autocov(%3u) = %18.15f [dB], i.e. autocovariance of SINR with tau = %3u\n',i,autocovGammaT1(i),i)
            if autocovGammaT1(i) < 10 * eps %%%% how it could enter this if condition if it's inside the while loop
                autocovGammaT1(i) = 0;
            end
            i = i + 1;
        end
        fprintf('\n')        
        dcurr = maxNumConsecutiveDropouts;
        y1 = 1e3*epsb;
        d = dcurr;
        while y1 > epsb && d < dmax
            n = 10^5; % % % Parameter of the solver getmvncdf by Botev (2017).
            % % [Botev] The normal law under linear restrictions -- simulation and estimation via minimax tilting
            xl = -Inf * ones(1,d);
            xu = (valGammaEpsp-mu) * ones(1,d);
            
            SIGMA = toeplitz([sigma^2, autocovGammaT1(1:d-1)]);
 
            est=getmvncdf(xl,xu,SIGMA,n);
            % % toc
            y0 = normcdf(valGammaEpsp,mu,sigma);
            y0 = y0^d;
            y1 = est.prob;
            ye = est.relErr;
            y2 = est.upbnd;
            fprintf('For d = %3u:\n', d)
            fprintf(' P_ind    = %17.15f\n', y0)
            fprintf(' P_approx = %17.15f\n', y1)
            fprintf(' P_upbnd  = %17.15f\n', y2)
            fprintf(' RelErr   = %17.15f\n', ye)
            fprintf(' DeltaP   = %17.15f\n', y1-y0)
            fprintf(' dP_abs   = %17.15f\n\n', y1-epsb)
            fprintf('The execution performed on ')
            fprintf(datestr(now))
            fprintf('\n\n')
            if y1-epsb > 0
                % d1 = d;
                % d = d1 + (dmax - d1)/2;
                d = d + 1;
            end
        end
        maxNumConsecutiveDropouts = d;
        
    end
   
    fprintf('       eta_A = %9.3f       i.e. expected value of PER\n', expValPER)
    fprintf('   sigma_A^2 = %9.3f       i.e. variance of PER\n', variancePER)
    if extNumDropouts == 1
        fprintf('  ell_B^star = %5u.          i.e. max num of consecutive dropouts; it depends on eps_b\n\n',...
            maxNumConsecutiveDropouts)
    end    
    % % % --- Auxiliary functions ----------------------------------------------------------------------------
    function cse0t = getcse0(t)
        cse0t = se0^2*exp(-1/2*(t/te0)^2);        
    end
    function cse1t = getcse1(t)
        cse1t = se1^2*exp(-1/2*(t/te1)^2); 
    end
    function cb0t = getcb0(t)
        cb0t = sb0^2*exp(-1/2*(V0*t/dc0)^2);
    end
    function cb1t = getcb1(t)
        cb1t = sb1^2*exp(-1/2*(V1*t/dc1)^2);
    end
    function M2t = getM2t(T)
        t = T + 1;
        sy0t = se0^2+sb0^2;
        sy1t = se1^2+sb1^2;
        cy0t = getcse0(t)+getcb0(t);
        cy1t=getcse1(t)+getcb1(t);
        Xt = exp(sy0t+cy0t)*D^2;
        Yt = 2*B*D*exp(1/2*sy1t+sy0t+cy0t);
        Zt = B^2*exp(sy1t+cy1t+sy0t+cy0t);
        M2t = Xt+Yt+Zt;
    end

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    % % % Z.I. Botev, "The normal law under linear restrictions:                                        % % % 
    % % %   simulation and estimation via minimax tilting,"                                             % % % 
    % % %       J. Roy. Statistical Soc.: Series B (Statistical Methodology), vol. 79, no. 1,           % % %  
    % % %           pp. 125–148, 2017.                                                                  % % % 
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
end