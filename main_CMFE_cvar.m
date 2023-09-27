% Assignment 1
% Portfolio optimization under CVaR constraint
% Geoffroy Hebert-Emond,  Vincent Lariviere, Benjamin Viau
% Deadline: 27/09/2023

% This script has 4 sections.

% Section 1 is where data is manually inputed by the used to give the
% specification of their problems.

% Section 2 construct the matrices and vector necessary to the estimation
% of the minimum variance portfolio and then perform this minimization
% using a function defined in the file "min_variance_portfolio.m". For more
% information, see the notes in this file.

% Section 3 estimate calculate the Value at Risk and the Conditonal Value
% at Risk attached to the results of the minimum variance portfolio.

% Section 4 estimates the minimum CVaR portfolio using the methodology from
% Tyrrell Rockafellar & Stanislav Uryasev in their paper "Optimization of
% conditional value-at-risk". It calculates the weight of each instrument
% in the optimal portfolio as well as the VaR and the CVaR associated with
% it at given sample size and risk levels. It also calculates the
% difference in VaR and CVaR between the minimum variance and the minimum
% CVaR approaches.

%% 1: Input data
%Covariance matrix of stocks
cov_matrix =   [0.00324625 0.00022983 0.00420395;
                0.00022983 0.00049937 0.00019247;
                0.00420395 0.00019247 0.00764097];

%Mean return of stocks
mean_returns = [0.0101110 0.0043532 0.0137058];

%Risk level
betas = [0.9, 0.95, 0.99];

%Sample size
q = [1000 3000 5000 10000 20000];

%Constraint on expected loss/return
R = 0.011;

%% 2: Approach 1: Minimum Variance

%Weight constraint
e_Aeq = ones(1, 3);
e_beq = 1;

%Return constraint
r_beq = -R;
A = -mean_returns;

%Combine constraints
Aeq = [e_Aeq; A];
beq = [e_beq; r_beq];


% Estimating Weight

weight_vector = quadprog(cov_matrix, [], [], [], Aeq, beq);


%% 3: Estimating VaR et CVaR
x = weight_vector;

portfolio_variance = x' * cov_matrix * x;
portfolio_stddev = sqrt(x' * cov_matrix * x);
portfolio_return = mean_returns * -x;

% Define the equality constraints, sum of weights equal 1
Aeq = [1, 1, 1];
beq = 1;

% Define bounds for asset weights and returns constraint
asset_num = length(mean_returns);
ub = ones(asset_num, 1);
lb = zeros(asset_num, 1);
b = -R;

% Calculate the min-variance optimal portfolio and extract values
cov_m = cov_matrix;
[opt_var, opt_returns] = min_variance_portfolio(A, b, Aeq, beq, lb, ub, cov_m);

% Define betas for VaR and CVaR calculations
betas = [0.9 0.95 0.99];
[VaR,CVaR] = VaR_CVaR(betas, opt_returns, opt_var);

min_var_results = table(betas',VaR,CVaR, 'VariableNames', {'betas', 'VaR', 'CVaR'});


%% 4: Approach 2: Minimum CVaR approach

%Instantiate operational arrays and variables
results = zeros(length(betas)*length(q),9);
counter = 1;

%Run the function for all beta and sample size
for i = betas
    for k = q
        counter
        [beta_risk,sample_size,x1,x2,x3,value_at_risk,c_value_at_risk,iters,time] = RiskCalculator(mean_returns,cov_matrix,i,k);
        results(counter,:) = [beta_risk,sample_size,x1,x2,x3,value_at_risk,c_value_at_risk,iters,time];
        counter = counter + 1;
    end
end

%Calculating VaR and CVaR diff
results = array2table(results, 'VariableNames', {'betas','sample_size','SP500','GovBond','SmallCap','value_at_risk','c_value_at_risk','iters','time'});
result_table = join(results, min_var_results, 'keys', 'betas');

result_table.('VaR diff. %') = ((result_table.('VaR') - result_table.('value_at_risk'))./result_table.('value_at_risk'))*100;
result_table.('CVaR diff. %') = ((result_table.('CVaR') - result_table.('c_value_at_risk'))./result_table.('c_value_at_risk'))*100;

% Create table
column_names = {'Beta', 'Sample Size', 'SP500', 'GovBond', 'SmallCap', 'VaR', 'VaR diff. (%)', 'CVaR', 'CVaR diff. (%)','Iterations', 'Time'};
table5 = table(result_table.('betas'), result_table.('sample_size'), result_table.('SP500'), result_table.('GovBond'), ...
    result_table.('SmallCap'), result_table.('value_at_risk'), result_table.('VaR diff. %'), result_table.('c_value_at_risk'), ...
    result_table.('CVaR diff. %'),result_table.('iters'),result_table.('time'), 'VariableNames', column_names);
            
disp(table5);
