

problem.objective = @root2d;
problem.x0 = 1;
problem.solver = 'fsolve';
problem.options = optimset(@fzero);
sol = fzero(problem);

  r=vpa(sol);
  r
   z=r(real(r)>0&imag(r)==0);
z
function F = root2d(x)

F(1) = (-x^N-x+2)/(x^3-x^(N+1)+x-1)-(1-p);
end
