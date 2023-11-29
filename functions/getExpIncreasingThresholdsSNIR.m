function SnirThresholds = getExpIncreasingThresholdsSNIR(nu,sigma,nChanLev,valGammaEpsp)
    
 p = quadgk(@getSteadyStateProbabilityIntervalSNIR,valGammaEpsp,Inf,'AbsTol',10*eps,'RelTol',100*eps);
    SnirThresholds = zeros(1,nChanLev-1);   % % SNIR's iterval thresholds
    p
    1-p
   % z=0;
%     syms x
%    equation=(-x^nChanLev-x+2)/(x^nChanLev-x^(nChanLev+1)+x-1)==(1-p);  
%    sol=solve(equation,x);
    z=my_exp_law(nChanLev,p);
%    problem.objective = @root2d;
%    problem.x0 = 1;
%    problem.solver = 'fsolve';
%    problem.options = optimset(@fzero);
%    sol = fzero(problem);
% 
%    r=vpa(sol);
%    z=r(real(r)>0&imag(r)==0);
   z
    %   invN = (1-p)/((2^(nChanLev-1))-1);
    %  invN=(1-p)/((3^(nChanLev-1))-2);
%         invN = 1/((z^(nChanLev))-1);
%         for index = 1 : nChanLev - 2
% %          SnirThresholds(index) = nu+sigma*sqrt(2)*erfinv(2*invN-1);
%            SnirThresholds(index) = norminv(invN,nu,sigma);        
%            invN= z*invN;
%         %  invN=(nChanLev-1)*invN;
%         
%         invN
% %          invN=invN+exp((index+1).*X);
%          
%         end
    %  i=1:N-2;
   
      for index = 1 : nChanLev - 2
      arg=((z.^index)-1)./(((z^nChanLev)-1).*(z-1));
     SnirThresholds(index) = norminv(arg,nu,sigma); 
      end
     
        
    SnirThresholds(nChanLev-1)= valGammaEpsp;
     
      function y = getSteadyStateProbabilityIntervalSNIR(x)
        y = normpdf(x,nu,sigma);
      end
%      function F = root2d(x)
% 
%         F(1) =(-x^nChanLev-x+2)/(x^nChanLev-x^(nChanLev+1)+x-1)-(1-p);
%      end

end
