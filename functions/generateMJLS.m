function [A_a,B_a,C_a,P_a,doubleState,received] = generateMJLS(A,B,Q,TPM,PER)
    % % % ----------------------------------------------------------------------------------------------------
    % % Generates a Markov jump linear system assosiated to the system controlled over the WHART channel
    % % % ----------------------------------------------------------------------------------------------------
    % %     Inputs:
    % %                    A:    system state matrix
    % %                    B:    system input matrix
    % %                    Q:    state weighting matrix in the optimal contriol problem
    % %                  TPM:    transition probability matrix (TPM) between the wireless channel's states
    % %                  PER:    packet error ratios (PER) in the wireless channel's states
    % % 
    % %     Outputs:
    % %         A_a:    sequence of system state matrices
    % %         B_a:    sequence of system input matrices
    % %        	C_a:    sequence of system otput matrices for the detectability test
    % %         P_a:    transition probability matrix of the MJLS, _a stands for augmented
    % % % ----------------------------------------------------------------------------------------------------
    numStates = size(TPM,2);
    doubleState = zeros(1,numStates);
    
    opModes = nnz(1 - PER) + nnz(PER);
    received = zeros(1,opModes);
    % % % ----------------------------------------------------------------------------------------------------
    index = 1;
    for theta = 1 : numStates
        if PER(theta) == 1
            received(index) = 0;
            index = index + 1;
        elseif PER(theta) == 0
            received(index) = 1;
            index = index + 1;
        else
            doubleState(theta) = 1;
            received(index) = 0;
            received(index+1) = 1;
            index = index + 2;
        end
    end
    % % % ----------------------------------------------------------------------------------------------------
    P_a = zeros(opModes); % % Initialization of the TPM of the MJLS, _A stands for Augmented
    ssDistribution = zeros(1,opModes);  % % steady state distribution
    offset_r = 0;           % % _r stands for row
    for theta_r = 1 : numStates
        if doubleState(theta_r) == 0
            % % %---------------------------------------------------------------------------------------------
            % % % Part for the computation of the steady_state_distribution:
            ssDistribution(theta_r+offset_r) = 1 / numStates;
            % % %---------------------------------------------------------------------------------------------
            offset_c = 0;   % % _c stands for column
            for theta_c = 1 : numStates
                if doubleState(theta_c) == 0       
                    P_a(theta_r+offset_r,theta_c+offset_c) = TPM(theta_r,theta_c);
                else
                    P_a(theta_r+offset_r,theta_c+offset_c) = PER(theta_c)*TPM(theta_r,theta_c);
                    offset_c = offset_c + 1;
                    P_a(theta_r+offset_r,theta_c+offset_c) = (1 - PER(theta_c))*TPM(theta_r,theta_c);
                end
            end
        else
            offset_c = 0;
            % % %---------------------------------------------------------------------------------------------
            % % % Part for the computation of the steady state distribution ssDistribution:
            ssDistribution(theta_r+offset_r) = PER(theta_r) / numStates;            
            % % %---------------------------------------------------------------------------------------------
            for theta_c = 1 : numStates
                if doubleState(theta_c) == 0       
                    P_a(theta_r+offset_r,theta_c+offset_c) = TPM(theta_r,theta_c);
                else
                    P_a(theta_r+offset_r,theta_c+offset_c) = PER(theta_c)*TPM(theta_r,theta_c);
                    offset_c = offset_c + 1;
                    P_a(theta_r+offset_r,theta_c+offset_c) = (1-PER(theta_c))*TPM(theta_r,theta_c);
                end
            end
            offset_r = offset_r + 1;
            % % %---------------------------------------------------------------------------------------------
            % % % Part for the computation of the steady state distribution:
            ssDistribution(theta_r+offset_r) = (1-PER(theta_r))/ numStates;  
            % % %---------------------------------------------------------------------------------------------
            P_a(theta_r+offset_r,:) = P_a(theta_r+offset_r - 1,:);
        end    
    end
    clear theta_r theta_c offset_r offset_c
    % % % ----------------------------------------------------------------------------------------------------
    A_a(:,:,opModes) = A;
    B_a(:,:,opModes) = B*received(opModes);
    C_a(:,:,opModes) = Q^(0.5);
    for i = 1 : opModes - 1
        A_a(:,:,i) = A;
        B_a(:,:,i) = B*received(i);
        C_a(:,:,i) = Q^(0.5);
    end 
end