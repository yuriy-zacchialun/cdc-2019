function [eta_M,var_M,nStarD] = getMarkovChannelLinkQualityMetrics(TPM,PEP,steadyProb,eps_b)
    % % % ----------------------------------------------------------------------------------------------------        
    fprintf('\n')
    eta_M = PEP*steadyProb';
    fprintf('       eta_M = %9.3f       i.e. expected value of PER over Markov channel\n', eta_M)
    var_M = PEP.^2*steadyProb' - eta_M^2;
    fprintf('   sigma_M^2 = %9.3f       i.e. variance of PER  over Markov channel\n', var_M)
    % nStarD = ceil(log(eps_b)/log(max(eig(TPM(1:end-1,1:end-1)))));
    nStarD = floor(log(eps_b)/log(max(eig(TPM(1:end-1,1:end-1)))));
    fprintf('      nStarD = %5u.          i.e. max num of consecutive dropouts; it depends on eps_b\n\n', ...
        nStarD)
    % % % ----------------------------------------------------------------------------------------------------
end