function  x_optimal = Resampling_MVO(mu, Q,targetRet, simulations, no_rep)
    n=size(Q,1);
    sample_data = cell(no_rep, 1);
    weights = cell(no_rep, 1);
    Variance = cell(no_rep, 1);
    means = zeros(no_rep, n);
    total_weights = zeros(n,1);  %initial the matrix for random generated data
    
    
    for i=1:no_rep
        random_ob = mvnrnd(mu, Q, simulations);
        sample_data{i} = random_ob;
        means(i, :) = geomean(1 + random_ob) - 1;
        Variance{i} = cov(random_ob);
        
        A = -means(i,:); 
        b = -targetRet;% orginally Ax>R, so -Ax<-R, therefore A and b both negative
  
        Aeq = ones(1,n);
        beq = 1;          % all variables sum up to one

        weights{i}=quadprog(Variance{i},[],A,b,Aeq,beq, [], []);
        total_weights=total_weights+weights{i}; % add up the weights
    end
    x_optimal=total_weights/no_rep;
end

