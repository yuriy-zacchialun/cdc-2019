function steadyStateProbabilitySNIR = getSteadyStateProbabilitySNIR(nu,sigma,nChanLev,SnirThreshold)    
    steadyStateProbabilitySNIR = zeros(1,nChanLev);   % % SNIR's iterval thresholds
    steadyStateProbabilitySNIR(1) = ...
        quadgk(@getSteadyStateProbabilityIntervalSNIR,-Inf,SnirThreshold(1),'AbsTol',10*eps,'RelTol',100*eps);
    for i = 2 : nChanLev - 1
        steadyStateProbabilitySNIR(i) = ...
            quadgk(@getSteadyStateProbabilityIntervalSNIR,SnirThreshold(i-1),SnirThreshold(i),...
                'AbsTol',10*eps,'RelTol',100*eps);
    end
    steadyStateProbabilitySNIR(nChanLev) = ...
        quadgk(@getSteadyStateProbabilityIntervalSNIR, SnirThreshold(nChanLev-1), Inf, ...
                'AbsTol',10*eps, 'RelTol',100*eps);
            %%SnirThreshold(nChanLev-1)= valGammaEpsp in the code of SnirThreshold
    if sum(steadyStateProbabilitySNIR) ~= 1
        % % % In case of numerical errors, we normalize the solution, rendering it stochastic:
        stochasticRow = steadyStateProbabilitySNIR/sum(steadyStateProbabilitySNIR);
        steadyStateProbabilitySNIR = stochasticRow;
    end   
    function y = getSteadyStateProbabilityIntervalSNIR(x)
        y = normpdf(x,nu,sigma);
    end
end