function [K_B, X_B] = getBernoullianController(n_x,A,B,Q,R,p_C,consoleMessageOff)
    % % % ----------------------------------------------------------------------------------------------------
    % % Finds the optimal state-feedback controller under TCP-like protocol and Bernoullian characterisation 
    % %     of the packet lossess, which is obtained from the unique stabilizing solution X_M of the
    % %     modified algebraic Riccati equation with characteristic parameter p_C. See
    % %     [Chen et al, 2016], and [Schenato et al, 2007] for additional details on the problem.
    % % % ----------------------------------------------------------------------------------------------------
    % %     Inputs:
    % %             n_x:    number of system state variables
    % %               A:    system state matrix
    % %               B:    system input matrix
    % %               Q:    state weight matrix in LQR problem
    % %               R:    control weight matrix in LQR problem
    % %             p_C:    probability of receiving a control packet (of Bernoullian i.i.d. random variable)
    % %             consoleMessageOff: 1 suppresses the trace of execution of the solver, 0 shows it       
    % %     Outputs:
    % %             K_B:    optimal gain matrix for TCP-like protocol and Bernoullian packet loss model
    % %             X_B:    infinite horizon solution of the associated modified algebraic Riccati equation 
    % % % ----------------------------------------------------------------------------------------------------
    
    % % % ----------------------------------------------------------------------------------------------------
    % % % Description of LMI system:
    % % % ----------------------------------------------------------------------------------------------------
    setlmis([])     % % % Initialize description of LMI system
    % % % ----------------------------------------------------------------------------------------------------
    % % % Matrix variables in LMI problem via lmivar(type,struct):
    % % % ----------------------------------------------------------------------------------------------------
    X = lmivar(1,[n_x 1]); % % % (symmetric, 1, of size n_x, full, 1, i.e., arbitrary symmetric matrix)
    % % % Specify term content of LMIs via lmiterm(termID,A,B,flag):
    lmiterm([-1 1 1 X],-1,1);
    lmiterm([-1 1 1 0],Q);
    lmiterm([-1 1 1 X],A',A);
    lmiterm([-1 1 2 X],p_C*A',B);
    lmiterm([-1 2 2 0],R);
    lmiterm([-1 2 2 X],B',B);
    lmiterm([-2 1 1 0],R);
    lmiterm([-2 1 1 X],B',B);
    lmiterm([-3 1 1 X],1,1);
    lmisys = getlmis;
    % % % ----------------------------------------------------------------------------------------------------
    % % % Solve the optimization problem via mincx, that minimizes linear objective under LMI constraints
    % % % ----------------------------------------------------------------------------------------------------
    % % % defcx helps to specify the objective function for mincx solver
    numVariables = decnbr(lmisys); % % % Total number of decision variables in system of LMIs
    c = zeros(1,numVariables);
for i = 1 : numVariables
    S = defcx(lmisys,i,X);
    c(i) = sum(diag(S));
end
    % % [tmin,xfeas] = feasp(lmisys);
    % fprintf('The infinite horizon optimal linear quadratic regulator is obtained by using a\n')
    % tstart = tic;             % % start a stopwatch timer to measure performance
    [~,xopt] = mincx(lmisys,-c,[1e-12,1e3,0,0,consoleMessageOff]); % % % Minimizes a linear objective under LMI constraints
    % telapsed = toc(tstart);     % % saves a stopwatch timer to measure performance
    % fprintf('getBernoullianController: Time spent by LMI solver: %8.2f [s]\n\n', telapsed) 
    X_B = dec2mat(lmisys,xopt,X);         % % _B stands for Bernoullian
    K_B = (B'*X_B*B + R)^(-1) * B'*X_B*A;   % % K_B is a Bernoullian optimal controller
    % % % ----------------------------------------------------------------------------------------------------
    % % [Chen et al, 2016] Chen, Michael Z.Q., and Zhang, Liangyin, and Su, Housheng, and Chen, Guanrong.
    % %     Stabilizing Solution and Parameter Dependence of Modified Algebraic Riccati Equation With 
    % %         Application to Discrete-Time Network Synchronization. IEEE Transactions on Automatic Control, 
    % %         Vol. 61, No. 1, January 2016.
    % % % ----------------------------------------------------------------------------------------------------
    % %     [Schenato et al, 2007] Schenato, Luca, and Sinopoli, Bruno, and Franceschetti, Massimo, and 
    % %         Poolla, Kameshwar, and Sastry, S. Shankar. 
    % %             Foundations of Control and Estimation Over Lossy Networks.
    % %                 Proceedings of the IEEE, Vol. 95, No. 1, January 2007. (pg. 171)
    % % % ----------------------------------------------------------------------------------------------------
end