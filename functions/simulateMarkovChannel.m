function [ack, theta] = simulateMarkovChannel(T, numStatesChannel, TPM_C, PER_C)
    % % % ----------------------------------------------------------------------------------------------------
    % % Simulate a possible behaviour of the Markov channel at the plant's sampling times 
    % % % ----------------------------------------------------------------------------------------------------
    % %     Inputs:
    % %                        T:   time horizon for simulation and plotting
    % %         numStatesChannel:   number of states in the Markov channel model of WirelessHART
    % %                    TPM_C:   Transition probability matrix (TPM) between Channel's states
    % %                    PER_C:   packet error ratios (PER) in the Channel's states
    % %     Outputs:
    % %               ack:  sequence of acknowledgements, indicates whether control packet was received
    % %             theta:  sequence of wireless channel states
    % % % ----------------------------------------------------------------------------------------------------
    ack = zeros(1,T);		 % % sequence of acknowledgements, indicates whether control packet was received
    theta = zeros(1,T+1);    % % sequence of wireless channel states
    % % % ----------------------------------------------------------------------------------------------------
    theta(1) = randi(numStatesChannel);
    for i = 2:T+1
        % % % ------------------------------------------------------------------------------------------------
        % %   Cumulative sum of a given row of TPM produces a cumulative distribution function (CDF).
        % % % ------------------------------------------------------------------------------------------------
        CDF = cumsum(TPM_C(theta(i-1),:));
        r = rand();
        for j = 1:numStatesChannel
            if r < CDF(j)
                theta(i) = j;
                break;
            end
        end
        r = rand();
        if r > PER_C(j)
        	ack(i-1) = 1;
        else
        	ack(i-1) = 0; 
        end       
    end
    % % % ----------------------------------------------------------------------------------------------------
end