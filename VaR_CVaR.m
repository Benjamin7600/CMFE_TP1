function [VaR, CVaR] = VaR_CVaR(betas, p_return, p_variance)
    %This script takes input in the form of an vector of probability level
    %for loss, a scalar of portfolio return and a scalar of portfolio
    %variance. It then calculates and output the CVaR and the VaR for each
    %probability level.

    % Example:
    % betas = [0.9 0.95]
    % p_return = 0.25
    % p_variance = 0.43
    % VaR_CVaR(betas, p_return, p_variance)
    
    % initialize VaR and CVaR variables
    VaR = zeros(length(betas), 1);
    CVaR = zeros(length(betas), 1);
        
    % Define portfolio std
    p_std = sqrt(p_variance);
        
    for x = 1:length(betas)
        beta = betas(x);
        n = norminv(1 - beta, 0, 1);
        
        % Calculate VaR for given beta
        VaR(x) = p_return - (n * p_std);
        
        %Calculate CVaR for given beta
        CVaR(x) = p_return - (1 / (1 - beta)) * -normpdf(n) * p_std;
    end
    risk_level = {'0.90'; '0.95'; '0.99'};

    % Create table for date
    betas_str = {'Risk Level (B)','VaR', 'CVaR'};
    T = table(risk_level, VaR, CVaR, 'VariableNames', betas_str);
        
    disp('Var and CVar obtained with min-var approach: ');
    disp(T);
end