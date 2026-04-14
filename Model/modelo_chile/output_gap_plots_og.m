% Input: M is an nxm matrix, where n is the number of observations (rows), and m is the number of variables (columns).
lambda = 1600;  % Smoothing parameter for quarterly data (use 1600 for quarterly data, 100 for annual data)

% Preallocate matrices for the output gap (cyclical component) and trend 
[n, m] = size(Ymat); 
v_output_gap = zeros(n, m);  % Cyclical component (output gap)
trend = zeros(n, m);       % Trend component

for j = 1:m
    % Apply HP filter to each column (variable)
    trend(:,j) = hpfilter(Ymat(:,j), lambda);
    v_output_gap(:,j) = Ymat(:,j) - hpfilter(Ymat(:,j), lambda); 
end

% Display the results for one of the variables
disp('Output gap (cyclical component) for the first variable:');
disp(v_output_gap(:, 1));

figure;
% Plot original series and trends
subplot(2, 1, 1); % Upper part for the original series and trends 
hold on;  
for j = 1:m
    plot(Ymat(1:20, j), 'DisplayName', ['Original Series - ' names{j}]);
    plot(trend(1:20, j), '--'); 
end 
hold off; 
legend off; 
title('Original Series and Trends for All Variables (First 20 Periods)'); 
xlabel('Time'); 
ylabel('Value');

% Plot output gaps
subplot(2, 1, 2); % Lower part for the output gaps 
hold on; 
for j = 1:m
    plot(v_output_gap(1:20, j), 'DisplayName', [names{j}]);
end 
hold off; 
legend; 
title('Output Gaps for All Variables (First 20 Periods)'); 
xlabel('Time'); 
ylabel('Output Gap');