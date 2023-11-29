function stabilizable = testStabilizability(A,B,P,consoleMessageOff)
% % % ----------------------------------------------------------------------------------------------------
% % Test of the stabilizability in PTI MJLs are based on [Costa, Marques and Fragoso, 2005], pg. 57--58.
% %     These conditions, extended here for all vertices, are necessary (not sufficient), i.e.
% %     if LMIs are unfeasible, the system is undetectable, otherwise an additional test on JSR is
% %     performed, to ensure the sufficiency and correctness.
% % % ----------------------------------------------------------------------------------------------------
[n_x,~,~] = size(A);        % % where n_x is the number of system state variables in a system state vector
[~,n_u,~] = size(B);        % % where n_u is the number of control input variables in a control vector
[~,opModes,V] = size(P);    % % where opModes is the number of the operational modes
% %       V is the number of vertices of covex polytope of TPMs
% % % ----------------------------------------------------------------------------------------------------
W1{opModes} = [];
W2{opModes} = [];
W3{opModes} = [];
% % % ----------------------------------------------------------------------------------------------------
% % % Specify the description of LMI system:
% % % ----------------------------------------------------------------------------------------------------
setlmis([])     % % % Initialize description of LMI system
% % % Specify matrix variables in LMI problem via lmivar(type,struct):
for i = 1 : opModes
    W1{i} = lmivar(1,[n_x 1]);  % % (symmetric, 1, of size n_x, full, 1, i.e., arbitrary symmetric matrix)
    W2{i} = lmivar(2,[n_x n_u]);% % (rectangular, 2, n_x-by-n_u matrix)
    W3{i} = lmivar(1,[n_u 1]);  % % (symmetric, 1, of size n_x, full, 1, i.e., arbitrary symmetric matrix)
end
% % % Specify term content of LMIs via lmiterm(termID,A,B,flag):
for i = 1 : opModes
    lmiterm([-i 1 1 W1{i}],1,1);
    lmiterm([-i 1 2 W2{i}],1,1);
    lmiterm([-i 2 2 W3{i}],1,1);
end
for i = 1 : opModes
    lmiterm([-i-opModes 1 1 W1{i}],1,1);
end
for l = 1 : V       % % for each vertex of the TPM
    for j = 1 : opModes
        for i = 1 : opModes
            lmiterm([j+(1+l)*opModes 1 1 W1{i}],P(i,j,l)*A(:,:,i),A(:,:,i)');
            lmiterm([j+(1+l)*opModes 1 1 W2{i}],P(i,j,l)*A(:,:,i),B(:,:,i)','s');
            lmiterm([j+(1+l)*opModes 1 1 W3{i}],P(i,j,l)*B(:,:,i),B(:,:,i)');
        end
        lmiterm([j+(1+l)*opModes 1 1 W1{j}],-1,1);
    end
end
% % % ----------------------------------------------------------------------------------------------------
lmisys = getlmis;
% % % ----------------------------------------------------------------------------------------------------
[~,xfeas] = feasp(lmisys,[0,0,0,0,consoleMessageOff]);

X{opModes} = zeros(n_x);
Y{opModes} = zeros(n_x,n_u);
for i = 1 : opModes
    X{i} = dec2mat(lmisys,xfeas,W1{i});
    Y{i} = dec2mat(lmisys,xfeas,W2{i});
end
% % % ----------------------------------------------------------------------------------------------------
K{opModes} = zeros(n_u,n_x);
for i = 1 : opModes
    K{i} = Y{i}'/X{i};
end
% % % ----------------------------------------------------------------------------------------------------
M10 = [];
for i = 1 : opModes
    M10 = blkdiag(M10, kron(A(:,:,i)+B(:,:,i)*K{i}, A(:,:,i)+B(:,:,i)*K{i}));
end
M0 = zeros(opModes*n_x^2, opModes*n_x^2, V);
for i = 1 : V
    M0(:,:,i) = kron(P(:,:,i),eye(n_x^2));
end
Lambda = zeros(opModes*n_x^2, opModes*n_x^2, V);
sr = zeros(1,V);
for i = 1 : V
    Lambda(:,:,i) = M0(:,:,i) * M10;
    sr(i) = max(abs(eig(Lambda(:,:,i))));
    disp(['testStabilizability: Spectral radius of Lambda_' num2str(i, 1) ' is ' num2str(sr(i), 9) '.'])
    disp(' ')
end
M = cell(1,V);
for i = 1 : V
    M{1,i} = Lambda(:,:,i);
end
% % % ----------------------------------------------------------------------------------------------------
if V > 1
    OPTIONS = jsrsettings ('logfile', 1, 'jsr.triang', 1);
    [BOUNDS, ~] = jsr(M, OPTIONS);
else
    BOUNDS(2) = max(abs(eig(Lambda(:,:,1))));
end
if BOUNDS(2) < 1
    stabilizable = 1;
else
    stabilizable = 0;
end

end
% % % ----------------------------------------------------------------------------------------------------
% % References:
% %     [Costa, Marques and Fragoso, 2005]
% %         Costa, Oswaldo Luiz do Valle and Marques, Ricardo Paulino and Fragoso, Marcelo Dutra.
% %             Discrete-time Markov jump linear systems. Springer-Verlag, London, 2005.
% % % ----------------------------------------------------------------------------------------------------
