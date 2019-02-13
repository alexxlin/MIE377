function  x_optimal = Robust_MVO(periodReturns, Q, confidence_level)



n=size(Q,1);

N = size(periodReturns,1);

mu=geomean(periodReturns+1)-1;
Theta=diag(diag(Q)./N);
Q_ = [Q zeros(n,1);zeros(1,n+1)];

kappa = sqrt(chi2inv(confidence_level,n));
lambda=50;




% Setup Gurobi Model
clear model;
% model.varnames = [tickers(1:n) 'aux.Var'];
model.Q=sparse(lambda*Q_);
model.obj = [-mu kappa];   % Zeros since no linear variables in obj. function
model.modelsense = 'min';

% linear constraints
model.A = sparse([ones(1,n),0]);
model.rhs=1;
model.sense='=';
model.lb=[-Inf(1,n) 0];

% Quadratic Constraints
model.quadcon(1).Qc = sparse([Theta zeros(n,1); zeros(1,n) -1]);
model.quadcon(1).q  = zeros(n+1,1);
model.quadcon(1).rhs = 0;
model.quadcon(1).sense = '=';

% Setup Gurobi Parameters
clear params;
params.TimeLimit = 100;
params.OutputFlag = 0;

% Solve robust MVO
rob_results = gurobi(model,params);
 x_optimal = rob_results.x(1:n);
end