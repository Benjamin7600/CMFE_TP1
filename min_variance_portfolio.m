function [var, returns] = min_variance_portfolio(A, b, Aeq, beq, lb, ub, cov_m)
    %This function takes a serie of constraint matrices and vectors (see
    %the documentation on linprog for the format) as well as a covariance
    %matrix and calculated the weights that minimize the variance of a
    %portfolio.
    %It outputs weights, the variance and the returns attached to the
    %optimal portfolio.

    % Example:
    % A = [5    1;
    %      3    3/4;
    %      4   -1;
    %     -1/4 -4;
    %     -2   -3;
    %     -1    2]
    % b = [5 1 6 8 9 3]
    % Aeq = [1 1/40]
    % beq = 1/2
    % lb = [0 0]
    % up = [4 6]
    % cov_m = [0.35 0.02;
    %          0.02 0.35]
    % min_variance_portfolio(A, b, Aeq, beq, lb, ub, cov_m)


    % Define objective function coeff's
    H = cov_m * 2;    
            
    % Define portfolio weights
    w = quadprog(H, [], [], [], [Aeq; A], [beq; b], lb, ub, []);
            
    % Calculate optimal portfolio variance and returns
    var = w' * cov_m * w;
    returns = A * w;
            
    % Create table for asset weights
    asset_name = {'SP500', 'GovBond', 'SmallCap'};
    T = table(w(1), w(2), w(3), 'VariableNames', asset_name);
            
    disp("The M-V portfolio variance is: ");
    disp(var);
    disp("The M-V portfolio expected return is: ");
    disp(returns);
    disp("The portfolio weights are:")
    disp(T);
end