function [expVal, vari] = getMomentsPERCoded(mu,sigma,lo,hi,ellF)
    
    expVal = quadgk(@giveMomentsCoded, lo, hi, 'AbsTol',10*eps, 'RelTol',100*eps); 
    %quadgk(fun,a,b,param1,val1,param2,val2,...) performs the integration with specified values of optional parameters
    % AbsTol: Absolute error tolerance. RelTol: Relative error tolerance.
    expSquared = quadgk(@giveMomentsSquareCoded, lo, hi, 'AbsTol',10*eps, 'RelTol',100*eps);
    
    vari = expSquared - expVal^2;
    
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% We use RS(32,26) with 8-bit symbols
%%% (see e.g. 
%%%    https://en.wikipedia.org/wiki/Reed–Solomon_error_correction
%%%    https://www.cs.cmu.edu/~guyb/realworld/reedsolomon/reed_solomon_codes.html )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function y = giveMomentsCoded(x)
        f = normpdf(x,mu,sigma);
        % % % PER:
           % return
        fprintf('consider the packet error rate (PER) to be r_p = %18.15e.\n',epsp)
        n = 32*ellF;
        k = 26*ellF;
          t = (n-k)/2;
         fprintf('For the RS(n,k,n-k+1)_q code with n = %u, k = %u, q = %u\n\n',n,k,2^ellF)
          syms x z
           f1 = 0;
           z= 4*exp(-10*10.^(x/10)) + (182*exp(-15*10.^(x/10)))/3 - (728*exp(-16*10.^(x/10)))/5 + ...
            (4004*exp(-18*10.^(x/10)))/15 + 429*exp(-(35*10.^(x/10))/2) - (56*exp(-(40*10.^(x/10))/3))/3 + ...
            (4004*exp(-(50*10.^(x/10))/3))/15 + (182*exp(-(55*10.^(x/10))/3))/3 - ...
            (8*exp(-(56*10.^(x/10))/3))/15 + exp(-(75*10.^(x/10))/4)/30 - (1144*exp(-(120*10.^(x/10))/7))/3 + ...
            4*exp(-(130*10.^(x/10))/7) - (1144*exp(-(160*10.^(x/10))/9))/3 - ...
            (728*exp(-(200*10.^(x/10))/11))/5 - (56*exp(-(240*10.^(x/10))/13))/3 ;
          for i = t+1:n
           f1 = nchoosek(vpa(n),i) * (z)^i * (1-z)^(n-i) + f1 ;   
          end
        % % % PDR:
        % % % g = 1 - Y; y = g.*f;
        y = f1.*f;    
    end
    function y = giveMomentsSquareCoded(x)
        f = normpdf(x,mu,sigma);
        % % % PER:
               fprintf('consider the packet error rate (PER) to be r_p = %18.15e.\n',epsp)
        n = 32*ellF;
        k = 26*ellF;
          t = (n-k)/2;
         fprintf('For the RS(n,k,n-k+1)_q code with n = %u, k = %u, q = %u\n\n',n,k,2^ellF)
          syms x z
           f1 = 0;
           z= 4*exp(-10*10.^(x/10)) + (182*exp(-15*10.^(x/10)))/3 - (728*exp(-16*10.^(x/10)))/5 + ...
            (4004*exp(-18*10.^(x/10)))/15 + 429*exp(-(35*10.^(x/10))/2) - (56*exp(-(40*10.^(x/10))/3))/3 + ...
            (4004*exp(-(50*10.^(x/10))/3))/15 + (182*exp(-(55*10.^(x/10))/3))/3 - ...
            (8*exp(-(56*10.^(x/10))/3))/15 + exp(-(75*10.^(x/10))/4)/30 - (1144*exp(-(120*10.^(x/10))/7))/3 + ...
            4*exp(-(130*10.^(x/10))/7) - (1144*exp(-(160*10.^(x/10))/9))/3 - ...
            (728*exp(-(200*10.^(x/10))/11))/5 - (56*exp(-(240*10.^(x/10))/13))/3 ;
          for i = t+1:n
           f1 = nchoosek(vpa(n),i) * (z)^i * (1-z)^(n-i) + f1 ;   
          end
        % % % PDR:
        % % % g = 1 - Y; y = g.^2.*f;
        y = f1.^2.*f;
    end
end