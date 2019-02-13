function  x_optimal = MVO(mu, Q, targetRet)
    
    % Use this function to construct your MVO portfolio subject to the
    % target return, with short-selling allowed. 
    %
    % You may use quadprog, Gurobi, or any other optimizer you are familiar
    % with. Just be sure to comment on your code to (briefly) explain your
    % procedure.

    % Find the total number of assets
    n = size(Q,1); 

    % *************** WRITE YOUR CODE HERE ***************
    %----------------------------------------------------------------------
    
    A = -mu'; 
    b = -targetRet;% orginally Ax>R, so -Ax<-R, therefore A and b both negative
   
    Aeq = ones(1,n);
    beq = 1;          % all variables sum up to one

    [X FVAL]=quadprog(Q,[],A,b,Aeq,beq);
    
    
    
    
    
    
    
    x_optimal = X;
    
    %----------------------------------------------------------------------
    
end