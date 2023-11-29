function [expVal, vari] = getMomentsPER(mu,sigma,lo,hi,ellF)
    
    expVal = quadgk(@giveMoments, lo, hi, 'AbsTol',10*eps, 'RelTol',100*eps); 
    %quadgk(fun,a,b,param1,val1,param2,val2,...) performs the integration with specified values of optional parameters
    % AbsTol: Absolute error tolerance. RelTol: Relative error tolerance.
    expSquared = quadgk(@giveMomentsSquare, lo, hi, 'AbsTol',10*eps, 'RelTol',100*eps);
    
    vari = expSquared - expVal^2;

    function y = giveMoments(x)
        f = normpdf(x,mu,sigma);
        % % % PER:
        Y = 1 - (4*exp(-10*10.^(x/10)) + (182*exp(-15*10.^(x/10)))/3 - (728*exp(-16*10.^(x/10)))/5 + ...
            (4004*exp(-18*10.^(x/10)))/15 + 429*exp(-(35*10.^(x/10))/2) - (56*exp(-(40*10.^(x/10))/3))/3 + ...
            (4004*exp(-(50*10.^(x/10))/3))/15 + (182*exp(-(55*10.^(x/10))/3))/3 - ...
            (8*exp(-(56*10.^(x/10))/3))/15 + exp(-(75*10.^(x/10))/4)/30 - (1144*exp(-(120*10.^(x/10))/7))/3 + ...
            4*exp(-(130*10.^(x/10))/7) - (1144*exp(-(160*10.^(x/10))/9))/3 - ...
            (728*exp(-(200*10.^(x/10))/11))/5 - (56*exp(-(240*10.^(x/10))/13))/3 - 1).^ellF;
        % % % PDR:
        % % % g = 1 - Y; y = g.*f;
        y = Y.*f;    
    end
    function y = giveMomentsSquare(x)
        f = normpdf(x,mu,sigma);
        % % % PER:
        Y = 1 - (4*exp(-10*10.^(x/10)) + (182*exp(-15*10.^(x/10)))/3 - (728*exp(-16*10.^(x/10)))/5 + ...
            (4004*exp(-18*10.^(x/10)))/15 + 429*exp(-(35*10.^(x/10))/2) - (56*exp(-(40*10.^(x/10))/3))/3 + ...
            (4004*exp(-(50*10.^(x/10))/3))/15 + (182*exp(-(55*10.^(x/10))/3))/3 - ...
            (8*exp(-(56*10.^(x/10))/3))/15 + exp(-(75*10.^(x/10))/4)/30 - (1144*exp(-(120*10.^(x/10))/7))/3 + ...
            4*exp(-(130*10.^(x/10))/7) - (1144*exp(-(160*10.^(x/10))/9))/3 - ...
            (728*exp(-(200*10.^(x/10))/11))/5 - (56*exp(-(240*10.^(x/10))/13))/3 - 1).^ellF;
        % % % PDR:
        % % % g = 1 - Y; y = g.^2.*f;
        y = Y.^2.*f;
    end
end