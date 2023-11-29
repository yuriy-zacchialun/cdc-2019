function [spectRadius, stabilizable] = decayConstant(A,B,P,K)
    [~,opModes,V] = size(P); 
    if iscell(K) == 0
        K1{opModes} = [];
        for i = 1 : opModes
            K1{i} = K;
        end
        K = K1;
        clear K1
    end
    M10 = [];
    for i = 1 : opModes
       M10 = blkdiag(M10, kron(A(:,:,i)-B(:,:,i)*K{i}, A(:,:,i)-B(:,:,i)*K{i}));
    end
    n_x = size(A,1);
    M0 = zeros(opModes*n_x^2, opModes*n_x^2, V);
    for i = 1 : V
        M0(:,:,i) = kron(P(:,:,i),eye(n_x^2));
    end    
    Lambda = zeros(opModes*n_x^2, opModes*n_x^2, V);
    sr = zeros(1,V);
    for i = 1 : V
        Lambda(:,:,i) = M0(:,:,i) * M10;
        sr(i) = max(abs(eig(Lambda(:,:,i))));
        disp(['Spectral radius ' num2str(i, 1) ' is ' num2str(sr(i), 9)])
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
   spectRadius = BOUNDS(2);
end