function [TPM] = getMarkovChannelSymbolRateLCR(...
            nu,sigma,autocov0der2,nChanLev,symFreq0,SnirThreshold,steadyStateProbabilitySNIR)
    % % % ----------------------------------------------------------------------------------------------------
    % % Generate Markov chain for a Gaussian with nonnull mean
    % % % ----------------------------------------------------------------------------------------------------
    % % Initialization:
    tPrev = zeros(1,nChanLev);              % % Transition probability of going in the previous state
    tNext = zeros(1,nChanLev);              % % Transition probability of going in the next state
    tSelf = zeros(1,nChanLev);              % % Transition probability of stating in the same state
    TPM = zeros(nChanLev,nChanLev);         % % Transition probability matrix
    % % % ----------------------------------------------------------------------------------------------------
    LCRnext = 1/(2*pi)*sqrt(autocov0der2/(sigma^2))*exp(-(SnirThreshold(1)-nu)^2/(2*sigma^2));
    RtiSym = symFreq0*steadyStateProbabilitySNIR;    
    % % % LCR analysis is valid only under the assumption that LCR << RtiSym. We verify this assumption as:
    orderMagnitude = 1e3;   % % % Threshold for the generation of the warning message
    if orderMagnitude*LCRnext > RtiSym(1)
        warning('The LCR analysis hypothesis are not satisfied!')
        fprintf('LCR/RtiSym = %5.3e\n\n',LCRnext/RtiSym(1))
    end    
    tPrev(1) = 0;
    tNext(1) = LCRnext/RtiSym(1);
    tSelf(1) = 1 - tNext(1);
    for j = 2 : nChanLev-1
        LCRcurr = LCRnext;
        LCRnext = 1/(2*pi)*sqrt(autocov0der2/(sigma^2))*exp(-(SnirThreshold(j)-nu)^2/(2*sigma^2));
        if orderMagnitude*LCRnext > RtiSym(j)
            warning('The LCR analysis hypothesis are not satisfied!')
            fprintf('LCR/RtiSym = %5.3e\n\n',LCRnext/RtiSym(j))
        end
        tPrev(j) = LCRcurr/RtiSym(j);
        tNext(j) = LCRnext/RtiSym(j);
        tSelf(j) = 1 - (tNext(j) + tPrev(j));
    end
    LCRcurr = LCRnext;
    tPrev(nChanLev) = LCRcurr/RtiSym(nChanLev);
    tNext(nChanLev) = 0;
    tSelf(nChanLev) = 1 - tPrev(nChanLev);    
    for i = 1 : nChanLev
        TPM(i,i) = tSelf(i);
    end
    for i = 2 : nChanLev
        TPM(i-1,i) = tNext(i-1);
    end
    for i = 1 : nChanLev-1
        TPM(i+1,i) = tPrev(i+1);
    end
    % % % If there are some substochastic rows in TPM, we render them stochastic:
    for i = 1 : size(TPM,2)
        TPM(i,:) = TPM(i,:)./sum(TPM(i,:));
    end
end