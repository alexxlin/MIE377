function  x_optimal = Diverse_MVO(mu, Q, targetRet,k)
n=size(Q,1);
Simularity_Mat=corrcov(Q); %Correlation Matrix
f = reshape(Simularity_Mat,[n*n,1]);
f = [f;zeros(n,1)]; % there are total of 20x20+20=420 assets

beq = [k;ones(n,1)];
first_row=[zeros(1,n*n) ones(1,20)];
Aeq_below=[];
for i=1:n
    temp=zeros(1,n*n+n);
    temp(((i-1)*n+1):((i-1)*n+n))=ones(1,n);
    Aeq_below=[Aeq_below;temp];
end  
Aeq=[first_row;Aeq_below]; %Setting Aeq and beq


I=eye(n*n);                                                                            
temp=-eye(n);
rightmatrix=repmat(temp,20,1);
A=[I,rightmatrix];
b=zeros(n*n,1);  % Setting A and b


 
lb=zeros(n*n+n,1);
ub=ones(n*n+n,1);
intcon=linspace(1,n*n+n, n*n+n); %lower & Upper bound, integer requirement


x=intlinprog(-f,intcon,A,b,Aeq,beq,lb,ub);


y=x(n*n+1:n*n+n,1); %The assets which have been chosen

A2 = -mu'; 
b2 = -targetRet;% orginally Ax>R, so -Ax<-R, therefore A and b both negative
lb2 = -inf*ones(1,n);
ub2 = inf*ones(1,n);  
for i=1:n
    if y(i)==0
        lb2(i)=0;
        ub2(i)=0;
    end
end
Aeq2 = ones(1,n);
beq2 = 1;          % all variables sum up to one
[X FVAL]=quadprog(Q,[],A2,b2,Aeq2,beq2,lb2,ub2);
x_optimal = X;  
end
    
    
    
    
    
    


       
        

