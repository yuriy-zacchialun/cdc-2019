function y=my_exp_law(N,p)

 syms x
 equation=(x.^(N-1)-1)./((x-1).*(x.^N)-1)-(1-p);  
 sol=solve(equation,x);
 
 sol_num=vpa(sol);
 y=max(sol_num(real(sol_num)>0&imag(sol_num)==0));
 y=double(y);

end