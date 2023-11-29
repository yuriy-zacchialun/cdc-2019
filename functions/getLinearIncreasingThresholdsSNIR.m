function SnirThresholds = getLinearIncreasingThresholdsSNIR(nu,sigma,nChanLev,valGammaEpsp)

 p = quadgk(@getSteadyStateProbabilityIntervalSNIR,valGammaEpsp,Inf,'AbsTol',10*eps,'RelTol',100*eps);
    SnirThresholds = zeros(1,nChanLev-1);   % % SNIR's iterval thresholds
    p
    1-p
   if mod(nChanLev,2) == 0
    invN_even = (1-p)/((nChanLev-1)*(nChanLev/2));
  
    for index = 1 : nChanLev - 2
        SnirThresholds(index) = nu+sigma*sqrt(2)*erfinv(2*invN_even-1);
        
        invN_even = invN_even+(index+1)*((1-p)/((nChanLev-1)*(nChanLev/2))); 
    end
    
   else
       invN_odd = (1-p)/((nChanLev)*((nChanLev-1)/2));

    for index = 1 : nChanLev - 2
        SnirThresholds(index) = nu+sigma*sqrt(2)*erfinv(2*invN_odd-1);
        invN_odd = invN_odd+(index+1)*((1-p)/((nChanLev)*((nChanLev-1)/2))); 
    end
       
       
   end
   
     SnirThresholds(nChanLev-1)= valGammaEpsp;
     
      function y = getSteadyStateProbabilityIntervalSNIR(x)
        y = normpdf(x,nu,sigma);
      end

end