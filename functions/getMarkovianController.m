function [K,J_av] = getMarkovianController(n_x,n_u,A,B,Q,R,Sigma_w,TPM,PER,consoleMessageOff) 
    %%%numStatesController,symFreq0,Ts,p_max,partitioningMethod,consoleMessageOff,ProbPacketError
    % % % ----------------------------------------------------------------------------------------------------
    % % Finds the optimal state-feedback controller under TCP-like protocol and Markovian characterisation 
    % %     of the packet lossess, which is obtained from the unique stabilizing solution X_{emme} of the
    % %     coupled algebraic Riccati equation with one time-step delayed mode observations. See [Matei, 2008] 
    % %     for details on the optimal LQR for MJLSs in presence of one time-step delayed mode observations.
    % % % ----------------------------------------------------------------------------------------------------
    % %     Inputs:
    % %                  n_x:    number of system state variables
    % %                  n_u:    number of control variables
    % %                    A:    system state matrix
    % %                    B:    system input matrix
    % %                    Q:    state weight matrix in LQR problem
    % %                    R:    control weight matrix in LQR problem
    % %              Sigma_w:    process noise covariance matrix
    % %    consoleMessageOff:    1 suppresses the trace of execution of the solver, 0 shows it 
    % % 
    % %     Outputs:
    % %         K:    optimal gain matrices for TCP-like protocol and Markovian packet loss model
    % %      J_av:    optimal cost of the Markovian controller, see [Matei, 2008]
    % % % ----------------------------------------------------------------------------------------------------
    numStatesController = length(TPM);
    doubleState = zeros(1,numStatesController);
    opModes = nnz(1 - PER) + nnz(PER);
    received = zeros(1,opModes);
    % % % ----------------------------------------------------------------------------------------------------
    index = 1;
    for theta = 1 : numStatesController
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
    TPM_A = zeros(opModes); % % Initialization of the TPM of the MJLS, _A stands for Augmented
    ssDistribution = zeros(1,opModes);  % % steady state distribution
    offset_r = 0;           % % _r stands for row
    for theta_r = 1 : numStatesController
        if doubleState(theta_r) == 0
            % % %---------------------------------------------------------------------------------------------
            % % % Part for the computation of the steady_state_distribution:
            ssDistribution(theta_r+offset_r) = 1 / numStatesController;
            % % %---------------------------------------------------------------------------------------------
            offset_c = 0;   % % _c stands for column
            for theta_c = 1 : numStatesController
                if doubleState(theta_c) == 0       
                    TPM_A(theta_r+offset_r,theta_c+offset_c) = TPM(theta_r,theta_c);
                else
                    TPM_A(theta_r+offset_r,theta_c+offset_c) = PER(theta_c)*TPM(theta_r,theta_c);
                    offset_c = offset_c + 1;
                    TPM_A(theta_r+offset_r,theta_c+offset_c) = (1 - PER(theta_c))*TPM(theta_r,theta_c);
                end
            end
        else
            offset_c = 0;
            % % %---------------------------------------------------------------------------------------------
            % % % Part for the computation of the steady state distribution ssDistribution:
            ssDistribution(theta_r+offset_r) = PER(theta_r) / numStatesController;            
            % % %---------------------------------------------------------------------------------------------
            for theta_c = 1 : numStatesController
                if doubleState(theta_c) == 0       
                    TPM_A(theta_r+offset_r,theta_c+offset_c) = TPM(theta_r,theta_c);
                else
                    TPM_A(theta_r+offset_r,theta_c+offset_c) = PER(theta_c)*TPM(theta_r,theta_c);
                    offset_c = offset_c + 1;
                    TPM_A(theta_r+offset_r,theta_c+offset_c) = (1-PER(theta_c))*TPM(theta_r,theta_c);
                end
            end
            offset_r = offset_r + 1;
            % % %---------------------------------------------------------------------------------------------
            % % % Part for the computation of the steady state distribution:
            ssDistribution(theta_r+offset_r) = (1-PER(theta_r))/ numStatesController;            
            % % %---------------------------------------------------------------------------------------------
            TPM_A(theta_r+offset_r,:) = TPM_A(theta_r+offset_r - 1,:);
        end    
    end
    clear theta_r theta_c offset_r offset_c
    % % % ----------------------------------------------------------------------------------------------------
    % % We find the Optimal Linear Quadratic Regulator for the augmented system, that is a Markov(ian) Jump 
    % %     Linear System (MJLS) with one time-step delayed mode observations. See [Matei, 2008] for details 
    % %         on the optimal LQR for MJLSs in the presence of one time-step delayed mode observations.    
    X{opModes} = [];
    % % % ----------------------------------------------------------------------------------------------------
    % % % Specify the description of LMI system:
    % % % ----------------------------------------------------------------------------------------------------
    setlmis([])     % % % Initialize description of LMI system
    % % % Specify matrix variables in LMI problem via lmivar(type,struct):
    for i = 1 : opModes
        X{i} = lmivar(1,[n_x 1]); % % % (symmetric, 1, of size n_x, full, 1, i.e., arbitrary symmetric matrix)
    end
    % % % Specify term content of LMIs via lmiterm(termID,A,B,flag):
for i = 1 : opModes
    lmiterm([-i 1 1 X{i}],-1,1);
    lmiterm([-i 1 1 0],Q);
    for j = 1 : opModes
        lmiterm([-i 1 1 X{j}],TPM_A(i,j)*A',A);
    end
    for j = 1 : opModes
        lmiterm([-i 1 2 X{j}],TPM_A(i,j)*A',B*received(j));
    end
    lmiterm([-i 2 2 0],R);
    for j = 1 : opModes
        lmiterm([-i 2 2 X{j}],TPM_A(i,j)*B',B*received(j));
    end
end
for i = 1 : opModes
    lmiterm([-i-opModes 1 1 0],R);
    for j = 1 : opModes
        lmiterm([-i-opModes 1 1 X{j}],TPM_A(i,j)*B',B*received(j));
    end
end
    lmisys = getlmis;
    % % % ----------------------------------------------------------------------------------------------------
    % % % Solve the optimization problem via mincx, that minimizes linear objective under LMI constraints
    % % % ----------------------------------------------------------------------------------------------------
    % % % defcx helps to specify the objective function for mincx solver
    numVariables = decnbr(lmisys); % % % Total number of decision variables in system of LMIs
    c = zeros(1,numVariables);
for i = 1 : numVariables
    S = defcx(lmisys,i,X{:});
    c(i) = sum(diag(S));
end
    % % % [tmin,xfeas] = feasp(lmisys,[0,0,0,0,consoleMessageOff]);
    % fprintf('The infinite horizon optimal linear quadratic regulator is obtained by using a\n')
    % tstart = tic;             % % start a stopwatch timer to measure performance
    [~,xopt] = mincx(lmisys,-c,[1e-16,1e3,0,0,consoleMessageOff]); % % % Minimizes a linear objective under LMI constraints
    X_opt{opModes} = [];
for i = 1 : opModes
    X_opt{i} = dec2mat(lmisys,xopt,X{i});
end
    % telapsed = toc(tstart);     % % saves a stopwatch timer to measure performance
    K_M{opModes} = [];
for i = 1 : opModes
    B_sum = zeros(n_u);
    C_sum = zeros(n_x,n_u);
    for j = 1 : opModes
        B_sum = TPM_A(i,j)*(B'*X_opt{j}*B*received(j) + R) + B_sum;
        C_sum = TPM_A(i,j)*(A'*X_opt{j}*B*received(j)) + C_sum;
    end
    K_M{i} = B_sum^(-1) * C_sum';
end
    % % % ----------------------------------------------------------------------------------------------------
    K{numStatesController}=[];
    offset_r = 0;
    for i = 1 : numStatesController
        if doubleState(i) > 0
            offset_r = offset_r + 1;        
        end
        K{i} = K_M{i+offset_r};
    end
    clear K_M lmisys i j c X S offset_r
    % % % ----------------------------------------------------------------------------------------------------
    % % % Optimal cost related to the considered controller according to [Matei et al, 2008] (a.k.a. J_av)
    % % % ----------------------------------------------------------------------------------------------------
    J_av = 0;
    for i = 1 : opModes
        J_av = J_av + ssDistribution(i) * trace(Sigma_w * X_opt{i});
    end
    % % %     fprintf('Summary for the considered optimal MJLQR with the number of states equal to\n')
    % % %     fprintf('\n%3u:\tp_R = %5.3f,\tJ_av  = %-5.3f, time spent by LMI solver: %8.2f [s]\n\n',...
    % % %         numStatesController,p_R,J_av,telapsed)  
    % % % ----------------------------------------------------------------------------------------------------
    % % References:
    % %     [Matei et al, 2008] Matei, Ion, and Martins, Nuno C., and Baras, John S.
    % %         Optimal Linear Quadratic Regulator for Markovian Jump Linear Systems, in the presence of 
    % %             one time-step delayed mode observations. Proceedings of the 17th World Congress of the
    % %                 International Federation of Automatic Control (IFAC), Seoul, Korea, July 6-11, 2008.
    % % % ----------------------------------------------------------------------------------------------------
end