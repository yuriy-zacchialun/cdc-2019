function TPM = getMarkovChannelSymbolRateAnalytic(nu,sigma,autocov1,threshold,steadyProb)            

    const2 = sigma^4-autocov1^2;
    const1 = 1/(2*pi*sqrt(sigma^4-autocov1^2));
    
    th = [-Inf, threshold, +Inf];
    nChanLev = length(steadyProb);
    TPM = zeros(nChanLev,nChanLev);
    
    for i = 1 : nChanLev
        for j = 1  : nChanLev
            P = integral2(@pdfSinrDouble, th(i), th(i+1), th(j), th(j+1), 'AbsTol',10*eps, 'RelTol',100*eps);
            TPM(i,j) = P/steadyProb(i);
        end
    end
    
    % % % If there are some nonstochastic rows in TPM, we render them stochastic:
    for i = 1 : size(TPM,2)
        TPM(i,:) = TPM(i,:)./sum(TPM(i,:));
    end
    
    function z = pdfSinrDouble(x,y)
        z = const1*exp(-0.5*(sigma^2*(x-nu).^2+sigma^2*(y-nu).^2-2*autocov1*(x-nu).*(y-nu))/const2); 
    end
end