function [stabile,rho] = testStabilityDelModeObs(A,B,K,P,nu)
    % % % ----------------------------------------------------------------------------------------------------
    % %     
    % % % ----------------------------------------------------------------------------------------------------
    N = size(P,1);
    n_x = size(A,1);
    deltaP1 = zeros(N^2,N);
    deltaPnu = zeros(N^2,N);
    bigA = zeros(n_x^2,N*n_x^2);
    bigBK = zeros(n_x^2,N*n_x^2);
    for i = 1 : N
    	deltaP1((i-1)*N+1:i*N,:) = diag(P(:,i));
    end
    for i = 1 : N
        bigA(:,(i-1)*n_x^2+1:i*n_x^2) = kron(conj(A),A);
    end
    for i = 1 : N
    	deltaPnu(:,i) = nu(i)*deltaP1(:,i);        
    end
    if iscell(K) == 0
        K1{N} = [];
        for i = 1 : N
            K1{i} = - K;
        end
        K = K1;
    else
        for i = 1 : N
            K{i} = - K{i};
        end
    end
    for i = 1 : N
        bigBK(:,(i-1)*n_x^2+1:i*n_x^2) = kron(conj(B*K{i}),B*K{i})+kron(conj(B*K{i}),A)+kron(conj(A),B*K{i});
    end
    Lambda = kron(deltaPnu,bigBK)+kron(deltaP1,bigA);
    rho = max(abs(eig(Lambda)));
    if rho < 1
        stabile = 1;
    else
        stabile = 0;
    end
end