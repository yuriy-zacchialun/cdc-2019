function PEP = getMarkovChannelPacketErrorProbability(nChanLev,mu,sigma,ellF,threshold,steadyProb,eps_p)  
    % % % ----------------------------------------------------------------------------------------------------
    PEP(nChanLev) = Inf;            % % % PEP = Packet Error Probability
    th = [-Inf, threshold, +Inf];
    
    for i = 1 : nChanLev
        [V, ~] = getMomentsPER(mu,sigma,th(i),th(i+1),ellF);
        if V < eps_p
            V = 0;  % % % The values under the threshold eps_p are considered pratically 0
        end
        PEP(i) = V/steadyProb(i);
    end    
    % % % ----------------------------------------------------------------------------------------------------
end