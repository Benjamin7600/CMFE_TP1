function [beta_risk,sample_size,x1,x2,x3,value_at_risk,c_value_at_risk,iters,time] = RiskCalculator(mean_return, cov_matrix, beta, sample_size)
    %This function takes as input a vector of mean return and a covariance
    %matrix from some asset to potentially include in the final portfolio
    %as well as a vector of risk level that we want to test at the level of
    %and a vector that defines the number of auxiliary variables.

    % Variables
    % mean_return: vector, shape: 1xn
    % cov_matrix: matrix, shape: nxn
    % beta: vector, shape: row vector of variable length
    % sample_size: vector, shape: row vector of variable length

    % Example:
    % mean_return = [0.015 0.26 0.34]
    % cov_matrix = [0.05  0.006 0.01
    %               0.006 0.24  0.014
    %               0.01  0.14  0.52]
    % beta = [0.95 0.99]
    % sample_size = [1000 5000]
    % RiskCalculator(mean_return, cov_matrix, beta, sample_size)

    tic
    q = sample_size;
    
    %Construct matrices
    %The matrices are constructed by creating submatrices that represent
    %the different components of the larger matrix and then fusing them
    %together. The columns vectors of the matrices are organized as follow:
    % [x1, x2, x3, alpha, u1, u2, ..., uk]
    % where xi represent the weight of asset i
    % alpha is the Value at Risk
    % ui are auxiliary variables
        
    %Objective function
    obj_return = 0*transpose(mean_return);
    obj_alpha  = 1;
    obj_u      = (1/(q*(1-beta)))*(ones(q,1));
        
    obj_function = [obj_return; obj_alpha; obj_u];
        
    %A
    A_y_simulation = mvnrnd(mean_return,cov_matrix,q);
    A_alpha = ones(size(zeros(q,1)));
    A_u = diag(ones(q,1));
        
    A = -horzcat(A_y_simulation,A_alpha,A_u);
        
    %b
    b = zeros(q,1);
        
    %Aeq
    Aeq_x_constraint      = ones(size(zeros(1,3)));
    Aeq_y_u_constraint    = zeros(1,q+1);
    Aeq_return_constraint = horzcat(-mean_return,zeros(1,q+1));
        
    sub_constraint = horzcat(Aeq_x_constraint,Aeq_y_u_constraint);
    Aeq = vertcat(sub_constraint,Aeq_return_constraint);
        
    %beq
    beq = [1;-0.011];
        
    %lower bound
    l_bound = vertcat(0,0,0,-inf,zeros(q,1));
        
    %higher bound
    h_bound = [];
        
    %options
    x0=optimset('MaxIter',10000000, 'algorithm', 'interior-point');
        
    %Estimating function
    [x, fval, output, iter] = linprog(obj_function,A,b,Aeq,beq,l_bound,h_bound,x0);
    
    runtime_func = toc;
    
    %Output
    beta_risk = beta;
    sample_size = sample_size;
    x1 = x(1);
    x2 = x(2);
    x3 = x(3);
    value_at_risk = x(4);
    c_value_at_risk = fval;
    iters = iter.iterations;
    time = runtime_func;

end