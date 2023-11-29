% % % ----------------------------------------------------------------------------------------------------
% % The main interface for the analysis of the TCP-like optimal networked state-feedback control over
% %     WirelessHART communication channel modelled as first order Markov chain
% % % ----------------------------------------------------------------------------------------------------
% % The following two commands should be executed always:
% % % ----------------------------------------------------------------------------------------------------
run('initialisation.m')     % % prepare the environment for execution
cd '/Users/yuriyzacchialun/Documents/MATLAB/cdc-2019'
run('config.m')             % % define all the relevant parameters listed in the configuration file

cd(folderData)
load('test.mat')

folderMain = '/Users/yuriyzacchialun/Documents/MATLAB/cdc-2019';
folderFigures =   [folderMain '/figures'];
folderData =      [folderMain '/data'];
folderFunctions = [folderMain '/functions'];
folderLog =       [folderMain '/log'];

% % % ----------------------------------------------------------------------------------------------------
% % Test Markov channel: stabilizability and detectability analysis
% % % ----------------------------------------------------------------------------------------------------    
if 1 - mu_M < p_max
    warning('The average probability of receiving a packet is too low!')
end
cd(folderFunctions)
[K_M,J_av] = getMarkovianController(n_x,n_u,A,B,Q,R,Sigma_w,mc_tp,mc_pe,~debugOn);
K_b = K_M;
for i = 1 : numStatesChannel
    K_b{i} = K_B;
end    
disp(' ')
[~,rho_B] =testStabilityDelModeObs(A,B,K_B,mc_tp,1-mc_pe);
fprintf('SpectRadius_Lambda_delayed_B   = %5.3f\n', rho_B)
[~,rho_M] =testStabilityDelModeObs(A,B,K_M,mc_tp,1-mc_pe);
fprintf('SpectRadius_Lambda_delayed_M   = %5.3f\n\n', rho_M)

numSimulations = 1e3; 

multiack = zeros(numSimulations,T);
multitheta = zeros(numSimulations,T+1);
multixm = zeros(n_x,T+1,numSimulations);
multiym = zeros(n_y,T+1,numSimulations);
multium = zeros(n_u,T,numSimulations);
multixb = zeros(n_x,T+1,numSimulations);
multiyb = zeros(n_y,T+1,numSimulations);
multiub = zeros(n_u,T,numSimulations);

clear generateNoise
w = generateNoise(n_x,Sigma_w,T);

fadingLength = 0;
maxFadingLength = 0;
for i = 1 : numSimulations
   [ack, theta] = simulateMarkovChannel(T, numStatesChannel, mc_tp, mc_pe);       
   fadingLength=max(accumarray(nonzeros((cumsum(ack)+1).*~ack),1));
   if maxFadingLength < fadingLength
       maxFadingLength = fadingLength;
   end
   [x_M, y_M, u_M] = simulateGilbertStateFeedbackControl(T,A,B,C,K_M,x_0,0*w,ack,theta);
   [x_B, y_B, u_B] = simulateGilbertStateFeedbackControl(T,A,B,C,K_b,x_0,0*w,ack,theta);
   multiack(i,:) = ack;
   multitheta(i,:) = theta;
   multixm(:,:,i) = x_M;
   multiym(:,:,i) = y_M;
   multium(:,:,i) = u_M;
   multixb(:,:,i) = x_B;
   multiyb(:,:,i) = y_B;
   multiub(:,:,i) = u_B;
end
maxym = zeros(n_y,T+1);
minym = zeros(n_y,T+1);
avgym = zeros(n_y,T+1);
maxyb = zeros(n_y,T+1);
minyb = zeros(n_y,T+1);
avgyb = zeros(n_y,T+1);
for i = 1 : n_y
    for j = 1 : T+1
        maxym(i,j) = max(multiym(i,j,:));
        minym(i,j) = min(multiym(i,j,:));
        avgym(i,j) = mean(multiym(i,j,:));            
        maxyb(i,j) = max(multiyb(i,j,:));
        minyb(i,j) = min(multiyb(i,j,:));
        avgyb(i,j) = mean(multiyb(i,j,:));
    end
end

fprintf('Max fading int considered in simulations lasts ')
fprintf('%u time steps.\n\n',maxFadingLength)

% % Performance analysis
% % % ----------------------------------------------------------------------------------------------------
J_B = trace(Sigma_w * X_B);
fprintf('\t\t\tJ_B   = %-5.3f\n',J_B)
fprintf('\t\t\tJ_M   = %-5.3f\n\n',J_av)
% % % ----------------------------------------------------------------------------------------------------
% % Comparison of the control gains
% % % ----------------------------------------------------------------------------------------------------
fprintf('    K_B   = [')
for i = 1 : n_x
    fprintf('%5.3f, ',K_B(i))
end
fprintf('\b\b]\n')
disp(' ')

for m = 1 : length(K_M)
    fprintf('    K_M(%02u) = [',m)
    for i = 1 : n_x
        fprintf('%5.3f, ',K_M{m}(i))
    end
    fprintf('\b\b]\n')
end

%% Plot the results of simulations
printOn = false;
cd(folderFigures)
plotSysOut(simHorizon,1/Rs,multiym,minym,maxym,avgym,plantType,...
    fullscreen,'outputMarkov',printOn)
plotSysOut(simHorizon,1/Rs,multiyb,minyb,maxyb,avgyb,plantType,...
    fullscreen,'outputBernoulli',printOn)

cd(folderMain)
fprintf('The execution terminated on ')
fprintf(datestr(now))
disp(' ')
diary OFF
