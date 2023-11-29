function [A,B,C,p_max,n_x,n_u,n_y] = getPendulum(M,m,b,I,g,l,Nbar,Ts,decimal,laggingAnalysis)
    % % % ----------------------------------------------------------------------------------------------------
    % % Provides the system matrices and upper bound on the critical probability of Bernulli channel for
    % % inverted pendulum on a cart, sampled at Ts, and computed from the system parameters presented in 
    % % [Franklin et al, 2009] G. F. Franklin, J. D. Powell, and A. Emami-Naeini. 
    % %     Feedback control of dynamic systems. Prentice Hall, 6th edition, 2009.
    % % It requires the Control System Toolbox to be installed.
    % % % ----------------------------------------------------------------------------------------------------
    % %     Inputs:
	% %               M:	mass of the cart                         [kg]
	% %               m:	mass of the pendulum                     [kg]
	% %               b:	coefficient of friction for cart         [N/m/sec]
	% %               I:	mass moment of inertia of the pendulum   [kg.m^2]
	% %               g:	gravitational acceleration                [m/s^2]
	% %               l:	length to pendulum center of mass        [m]
    % %              Ts:	sampling time                            [s]
    % %         decimal:    number of decimal digits of round(), that it is useful to present the results
    % % laggingAnalysis:    1 to perform analysis of the lagging effect associated with a 
    % %                     zero-order hold discretisation, 0 otherwise.
    % %     Outputs:
    % %                 A:  system state matrix
    % %                 B:  system input matrix
    % %                 C:  system output matrix
    % %             p_max:  upper bound on the critical probability of Bernulli channel
    % %               n_x:  number of system state variables
    % %               n_u:  number of control variables
    % %               n_y:  number of output variables
    % % % ----------------------------------------------------------------------------------------------------
    % % % Other quantities involved:
    % % %        force applied to the cart
    % % %        cart position coordinate
    % % %        pendulum angle from vertical
    % % % ----------------------------------------------------------------------------------------------------
    % % % The state variables of pendulum are 
    % % %      the cart’s position coordinate and pendulum’s angle from vertical,
    % % %      together with respective first derivatives.
    % % % ----------------------------------------------------------------------------------------------------
       
    % % % ----------------------------------------------------------------------------------------------------
    % % The state-space model of the inverted pendulum (on a cart).
    % % % ----------------------------------------------------------------------------------------------------
    % % The (linearized) state-space model is 
    % %     x_{k+1} = A * x_k + B * u_k
    % %         y_k = C * x_k
    % % where
    % %     x_k is a system state vector of n_x variables
    % %     u_k is a control vector of n_u variables
    % %     y_k is an output vector of n_y variables
    % % % ----------------------------------------------------------------------------------------------------
    
    % % % ----------------------------------------------------------------------------------------------------
    % % Continuos time state-space model:
    % % % ----------------------------------------------------------------------------------------------------
    p = I*(M+m)+M*m*l^2;                        % % denominator for the A and B matrices
    % % % ----------------------------------------------------------------------------------------------------
    cA = [0      1              0           0;
          0 -(I+m*l^2)*b/p  (m^2*g*l^2)/p   0;
          0      0              0           1;
          0 -(m*l*b)/p       m*g*l*(M+m)/p  0];
    % % % ----------------------------------------------------------------------------------------------------
    cB = [     0;
          (I+m*l^2)/p;
               0;
             m*l/p];
    % % % ----------------------------------------------------------------------------------------------------
    cC = [1 0 0 0;
          0 0 1 0];
    % % % ----------------------------------------------------------------------------------------------------
    sysc = ss(cA,cB*Nbar,cC,[]);
    % % % ----------------------------------------------------------------------------------------------------
    sys = c2d(sysc,Ts);
    % % % ----------------------------------------------------------------------------------------------------
    if laggingAnalysis == 1
    	impulse(sysc,3)         % % plots continuous output response
        hold on
        [x,t]=impulse(sys,3);
        plot(t,x,'ro')          % % plots discrete output response
        hold off
    end
    % % % ----------------------------------------------------------------------------------------------------
    A = round(sys.A,decimal);
    B = round(sys.B,decimal);
    C = round(sys.C,decimal);
    % % % ----------------------------------------------------------------------------------------------------
    % % disp(eig(A))
    dCo = ctrb(A,B);
    unco = length(A) - rank(dCo);
    if unco == 0
        fprintf('The described plant sampled at Ts = %5.3f is a controllable system.\n\n',Ts)
    else
        error('      The described plant sampled at Ts = %5.3f is uncontrollable!\n\n',Ts)       
    end
    % % % ----------------------------------------------------------------------------------------------------
    % clear b cA cB cC dCo g I l m M p sys sysc unco ell
    % % % ----------------------------------------------------------------------------------------------------
    fprintf('    A = [%7.3f, %7.3f, %7.3f, %7.3f;\n', A(1,1), A(1,2), A(1,3), A(1,4))
    fprintf('         %7.3f, %7.3f, %7.3f, %7.3f;\n', A(2,1), A(2,2), A(2,3), A(2,4))
    fprintf('         %7.3f, %7.3f, %7.3f, %7.3f;\n', A(3,1), A(3,2), A(3,3), A(3,4))
    fprintf('         %7.3f, %7.3f, %7.3f, %7.3f];\n\n',A(4,1), A(4,2), A(4,3), A(4,4))
    % % % ----------------------------------------------------------------------------------------------------
    fprintf('    B = [%7.3f;\n',    B(1,1))
    fprintf('         %7.3f;\n',    B(2,1))
    fprintf('         %7.3f;\n',    B(3,1))
    fprintf('         %7.3f];\n\n', B(4,1))
    % % % ----------------------------------------------------------------------------------------------------
    % % Find the limit probability of packet received for the Bernoulli representation of the chennel
    % %     following the procedure described in [Schenato et al, 2007].
    % % % ----------------------------------------------------------------------------------------------------
    eigA = eig(A);                  % % The eigenvalues of A
    p_min = 1 - 1/(max(eigA)^2);    % % The lower bound of the critical probability for Bernulli channel
    L = eigA > 1;                   % % All the unstable eigenvalues
    L = find(L);                    % % The indises of the unstable eigenvalues
    P = 1;
    for i = 1 : length(L)
        P = P * eigA(L(i))^2;
    end
    p_max = 1 - 1/P;            % % The upper bound of the critical probability for Bernulli channel
    fprintf('The lower bound p_min of the critical probability for Bernulli channel is %5.3f.\n', p_min)
    fprintf('The upper bound p_max of the critical probability for Bernulli channel is %5.3f.\n\n', p_max)
    % clear eigA p_min P L
    n_x = length(A);                % % number of system state variables
    n_u = size(B);
    n_u = n_u(2);                   % % number of control variables
    n_y = size(C);
    n_y = n_y(1);                   % % number of output variables   
end