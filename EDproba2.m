function p = EDproba2(v,s,k)
format shortG;
Xinit = zeros(4,1); %initial points for the fmincon function

options = optimoptions('fmincon',...
    'Algorithm','sqp-legacy','Display','iter','ConstraintTolerance',1e-30);

[xopt,fval,exitflag,output,lambda,grad,hessian] =...
    fmincon(@fct,Xinit,[],[],[],[],[],[],[],options);

function [f] = fct(Xinit)

N = 4;  % select number of variables

n=40; %order of Gauss-Legendre quadrature
[z,w]=legzo(n, -6,6); %function that returns weights and points of G-L quadrature

Mean = 0;
Variance = 1;
Skewness= s;
Kurtosis=  k;

A = z;
B= zeros(n,4);

B(:,1)= A';
B(:,2)= A'.^2;
B(:,3)= A'.^3;
B(:,4)= A'.^4;

 
for i = 1:n
C(i,1)= Mean;
C(i,2)= Variance;
C(i,3)= Skewness;
C(i,4)= Kurtosis;

end

D = B - C;

u = (D*Xinit);

f = sum((w)*exp(u)); %the function to minimize

end

%Compute the probability of v
p = exp((xopt(1)*(v-Mean) + xopt(2)*((v.^2) - Variance) + xopt(3)*((v.^3) - Skewness)...
    +xopt(4)*((v.^4) - Kurtosis))) ...
  /fval;


end