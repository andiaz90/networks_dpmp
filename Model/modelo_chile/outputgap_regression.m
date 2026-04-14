clc
% Example data: 
% Y is the dependent variable matrix (n x m) % X is the independent variable matrix (n x p) 
[n, m] = size(Ymat);  % n is the number of observations, m is the number of dependent variables 
[~, p] = size(pmat);  % p is the number of independent variables
lambda=1600;
% Preallocate matrix for storing regression coefficients 
beta = zeros(2, m+1);  % p + 1 because we include an intercept term

v_output_gap = zeros(n, m);  % Cyclical component (output gap)
pi_mat = diff((pmat));       % Trend component

for j = 1:m
    % Apply HP filter to each column (variable)
    v_output_gap(:,j) = Ymat(:,j)-hpfilter(Ymat(:,j), lambda); 
end

agg_y_gap=Yagg'-hpfilter(Yagg, lambda);
v_output_gap=[agg_y_gap  v_output_gap];
X=v_output_gap(1:end-1,:);
Y=[diff(pi)' pi_mat];
names2=[{'all'}; names];

% Loop through each column of Y and perform OLS regression 
    for j = 1:m+1
    % Add a column of ones to X to account for the intercept
    X_with_intercept = [ones(n-1, 1), X(:,j)];  % X_with_intercept is now (n x (p+1))
    
    % Perform OLS regression:
    % beta = (X'X)^(-1) * X'Y
    beta(:, j) = (X_with_intercept' * X_with_intercept) \ (X_with_intercept' * Y(:, j)); 
    end

% Create a table with the OLS coefficients var_names = strcat('Y', string(1:m));  % Variable names for the dependent variables in Y coef_names = ['Intercept'; strcat('X', string(1:p))];  % Names for the intercept and independent variables

% Convert the results into a table
coef_table = array2table(beta, 'VariableNames', names2, 'RowNames', {'alpha', 'beta'});

% Display the table
disp('OLS Regression Coefficients:');
disp(coef_table);

