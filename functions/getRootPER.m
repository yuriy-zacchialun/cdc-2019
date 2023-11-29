function valGammaEpsp = getRootPER(ellF,epsp)

    problem.objective = @(x)(getPER(x,ellF)-epsp)*1e6;
    % problem.objective = @(x)getPER(x,ellF)-epsp;
    problem.x0 = 1;
    problem.solver = 'fzero'; % a required part of the structure
    problem.options = optimset(@fzero); % default options
    % % % problem.options = optimset('Display','iter'); % show iterations
    %options = optimset(optimfun) creates an options structure options with all parameter names and ...
    %default values relevant to the optimization function optimfun.

    valGammaEpsp = fzero(problem); % % % Alternative: implement the Newton's method, with a step size control

    function y = getPER(x,ellF)
        y = 1 - (4*exp(-10*10.^(x/10)) + (182*exp(-15*10.^(x/10)))/3 - (728*exp(-16*10.^(x/10)))/5 + ...
            (4004*exp(-18*10.^(x/10)))/15 + 429*exp(-(35*10.^(x/10))/2) - (56*exp(-(40*10.^(x/10))/3))/3 + ...
            (4004*exp(-(50*10.^(x/10))/3))/15 + (182*exp(-(55*10.^(x/10))/3))/3 - ...
            (8*exp(-(56*10.^(x/10))/3))/15 + exp(-(75*10.^(x/10))/4)/30 - (1144*exp(-(120*10.^(x/10))/7))/3 + ...
            4*exp(-(130*10.^(x/10))/7) - (1144*exp(-(160*10.^(x/10))/9))/3 - ...
            (728*exp(-(200*10.^(x/10))/11))/5 - (56*exp(-(240*10.^(x/10))/13))/3 - 1).^ellF;
    end
end