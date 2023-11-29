function [x, y, u] = ...
    simulateGilbertStateFeedbackControl(T,A,B,C,K,x_0,w,ack,theta)
    % % % ----------------------------------------------------------------------------------------------------
    % %   Compute the state and output vectors for the time hotizon T applying Markovian control strategy
    % % % ----------------------------------------------------------------------------------------------------
    % %     Inputs:
    % %                   T:    time horizon for simulation and plotting
    % %                   A:    state matrix
    % %                   B:    input matrix
    % %                   C:    output matrix
    % %                   K:    optimal gain matrix for TCP-like protocol and Markovian packet loss model
    % %                 n_x:    number of system state variables
    % %                 n_y:    number of output variables
    % %                 x_0:    initial state
    % %                   w:    process noise
    % %                 ack:    sequence of acknowledgements, indicates whether control packet was received
    % %               theta:    sequence of wireless channel states
    % % numStatesController:    number of the operational modes of the controller
    % %                xt_R:    level crossing thresholds (xt) for the regulator states    
    % %     Outputs:
    % %                   x:    system state vector of n_x variables, T+1 time steps
    % %                   u:    system control inputs
    % % % ----------------------------------------------------------------------------------------------------
    
    % % % ----------------------------------------------------------------------------------------------------    
    % % The state-space model is 
    % %     x_{k+1} = A * x_k + ack_k * B_{M_k} * u_k + w_k
    % % where
    % %     x_k is a system state vector of n_x variables
    % %     u_k is a control vector of n_u variables, u_k = - K_{M_{k-1}} * x_k
    % %     w_k is a process noise with zero mean and covariance matrix Sigma_w
    % %   ack_k is a binary variable indicating whether control packet was received in time step k
    % % % ----------------------------------------------------------------------------------------------------
    % % The quadratic cost function for the system with multiple modes of operations in Markov Jump setting is
    % %     J_{av}(u) = limsup_{N \to Inf} ( 1/N * ...
    % %                     * E(sum_{k=0}^{N-1}((x_k' * Q * x_k + u_k' * R * u_k) + x_N' * Q * x_N))
    % %         where E() is the expextation, N < Inf, R and Q are positive definite matrices
    % % % ----------------------------------------------------------------------------------------------------
    [n_x,~] = size(A);        % % where n_x is the number of system state variables in a system state vector
    [~,n_u] = size(B);        % % where n_u is the number of system control inputs in a system contol vector
    [n_y,~] = size(C);        % % where n_y is the number of system output variables
    % % % ----------------------------------------------------------------------------------------------------
    % % Markov controller behaviour for the duration of the simulation
    % % % ----------------------------------------------------------------------------------------------------
    x = zeros(n_x,T+1);           % % state vector of n_x variables, T+1 time steps
    y = zeros(n_y,T+1);             % % output vector of n_y variables, T+1 time steps
    u = zeros(n_u,T);             % % control vector of n_u variables, T+1 time steps
    
    x(:,1) = x_0;
    x(:,2) = ( A - ack(1) * B * K{theta(1)} ) * x(:,1) + w(:,1); % % Under assumption of the availability 
                                      % % of the channel state information (CSI) emme(-1)=emme(0), as in 
                                      % % [Matei et al, 2008], see MarkovController.m for the reference
    u(:,1) = - K{theta(1)} * x(:,1);
    for l = 2:T
    	x(:,l+1) = ( A - ack(l) * B * K{theta(l-1)} ) * x(:,l) + w(:,l);
        u(:,l) = - K{theta(l-1)} * x(:,l);
    end
    for l = 1 : T+1
        y(:,l) = C * x(:,l);
    end
end