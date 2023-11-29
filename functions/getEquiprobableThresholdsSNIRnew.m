function SnirThresholds = getEquiprobableThresholdsSNIRnew(nu,sigma,nChanLev,valGammaEpsp)
%%%%% This function is to divide the bad state of Gilbert channel using
%%%%% equipropable partitioning method.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %p is the steadystate probability of the state where the PER=0 so the
    %performance is good
    p = quadgk(@getSteadyStateProbabilityIntervalSNIR,valGammaEpsp,Inf,'AbsTol',10*eps,'RelTol',100*eps);
    SnirThresholds = zeros(1,nChanLev-1);   % % SNIR's iterval thresholds
   
    invN = (1-p)/(nChanLev-1);

    for index = 1 : nChanLev - 2
        SnirThresholds(index) = nu+sigma*sqrt(2)*erfinv(2*invN-1);
        invN = invN+((1-p)/(nChanLev-1)); 
    end
     SnirThresholds(nChanLev-1)= valGammaEpsp;
     
      function y = getSteadyStateProbabilityIntervalSNIR(x)
        y = normpdf(x,nu,sigma);
      end
end