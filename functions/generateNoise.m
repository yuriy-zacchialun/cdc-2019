function w = generateNoise(n_x,Sigma_w,T)
    % % % ----------------------------------------------------------------------------------------------------
    % % Generates the process noise vector
    % %     Inputs:
    % %         n_x:        number of system state variables
    % %         Sigma_w:    process noise covariance matrix
    % %         T:          time horizon for simulation and plotting
    % % % ----------------------------------------------------------------------------------------------------
    mu = zeros(n_x,1);              % % n_x-dimensional 0-mean vector for the noise
    % % % ----------------------------------------------------------------------------------------------------
    w = mvnrnd(mu, Sigma_w, T+1);   % % noise vector of n_x variables, T+1 time steps
    w = w';
    % % % ----------------------------------------------------------------------------------------------------
end