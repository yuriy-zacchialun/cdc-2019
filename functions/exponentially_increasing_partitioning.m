% NO TRUNCATION
N=15;

i=1:N-1;

arg=((2.^i)-1)./((2^N)-1);
gamma=norminv(arg);


steady=zeros(1,N);

for j=2:N-1

    
steady(j)=normcdf(gamma(j))-normcdf(gamma(j-1));

end

steady(1)=normcdf(gamma(1));
steady(N)=1-sum(steady(1:N-1));

figure
stem(steady)

% TRUNCATION
N=10;
p_trunc=0.07;

syms y
p = ((y.^(N-1))-1)./((y-1).*((y.^N)-1))-(1-p_trunc);
R=solve(p,y);
Rnum=vpa(R);

steady=zeros(1,N);
i=1:N-2;


x_=Rnum(real(Rnum)>0&imag(Rnum)==0);
x=double(x_);
arg=((x.^i)-1)./(((x^N)-1).*(x-1));

mu=-5.235;
sigma=2.564;

gamma=norminv(arg,mu,sigma);


for j=2:N-2
   
steady(j)=normcdf(gamma(j),mu,sigma)-normcdf(gamma(j-1),mu,sigma);

end

steady(1)=normcdf(gamma(1),mu,sigma);
steady(N-1)=1-sum(steady(1:N-2))-p_trunc;
steady(N)=p_trunc;

figure
stem(steady)

