function [x_optimal, Ret] = cvar(mu, Q, targetRet,scenarios, S_0,confidence, dt)
n=size(Q,1);
Rho=corrcov(Q);
L=chol(Rho,'lower');
variance=diag(Q);
std_dev=sqrt(variance);
Paths=scenarios;
S=cell(1,n);
T=26 % six months
for i=1:n
    S{1,i}=[S_0(i)*ones(1,Paths)];
end

for i=1:Paths
    for j=1:T/dt
        Big_epsilon=L*randn(n,1);
        for q=1:n
            S{1,q}(j+1,i)=S{1,q}(j,i)*exp((mu(q)-1/2*variance(q))*dt+std_dev(q)*sqrt(dt)*Big_epsilon(q));
        end
    end
end



Ret=zeros(scenarios,20);
for i=1:Paths
    for j=1:n
        temp = T/dt+1;
        Ret(i,j)=S{1,j}(temp,i)/S{1,j}(1,i)-1;
    end
end
        
f=[1, 1/((1-confidence)*Paths)*ones(1,Paths), zeros(1,20)];
Aeq=[zeros(1,scenarios+1),ones(1,20)];

beq=1;
A=[-ones(scenarios,1), -eye(scenarios), -Ret];
b=[zeros(scenarios,1)];
lb=[-Inf;zeros(scenarios+20,1)];
ub=Inf*ones(scenarios+21,1);
[x,fval] = linprog(f,A,b,Aeq,beq,lb,ub);
x_optimal= x(scenarios+2:scenarios+21,1);
end
