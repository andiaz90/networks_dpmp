% Check I-O matrix orientation and values
clear; clc;

% Read the data
betaio = readmatrix('IO_2021_chile.csv');
betaio = betaio(1:12,1:12);

% Normalize by columns (each purchasing sector sums to 1)
betax = betaio ./ sum(betaio,1);

% Transpose so modbeta(i,j) = share that sector i purchases from sector j
modbeta = betax';

fprintf('=== I-O MATRIX STRUCTURE CHECK ===\n\n');

% Show dimensions
fprintf('Matrix dimensions: %dx%d\n\n', size(modbeta,1), size(modbeta,2));

% Check column sums (should be 1 after normalization and transpose)
col_sums = sum(modbeta,1);
fprintf('Column sums (each column = inputs purchased by one sector):\n');
for j = 1:12
    fprintf('  Sector %2d: %.6f\n', j, col_sums(j));
end
fprintf('\n');

% Check specific values
fprintf('=== KEY VALUES ===\n');
fprintf('modbeta(1,1) = %.6f (Sector 1 purchases from Sector 1)\n', modbeta(1,1));
fprintf('modbeta(1,3) = %.6f (Sector 1 purchases from Sector 3)\n', modbeta(1,3));
fprintf('modbeta(3,1) = %.6f (Sector 3 purchases from Sector 1)\n', modbeta(3,1));
fprintf('\n');

% Display first row (what sector 1 purchases from all sectors)
fprintf('First ROW of modbeta (Sector 1 as PURCHASER):\n');
fprintf('Shows what Sector 1 (Agriculture) buys from each sector:\n');
disp(modbeta(1,:));

% Display first column (what all sectors purchase from sector 1)
fprintf('First COLUMN of modbeta (Sector 1 as SUPPLIER):\n');
fprintf('Shows what each sector buys from Sector 1 (Agriculture):\n');
disp(modbeta(:,1)');

fprintf('\n=== INTERPRETATION ===\n');
fprintf('In the model: PM_i = sum_j beta(i,j) * PH_j\n');
fprintf('This means: beta(i,j) = share that sector i purchases from sector j\n');
fprintf('So modbeta(i,j) should have i=purchaser, j=supplier\n\n');

fprintf('In the TABLE: Element (i,j) shows "product i in sector j''s inputs"\n');
fprintf('This means: Table(i,j) has i=supplier, j=purchaser\n');
fprintf('So the table should show TRANSPOSE of modbeta!\n\n');

% What the table should show
table_matrix = modbeta';
fprintf('=== TABLE VALUES (should be transpose of modbeta) ===\n');
fprintf('First 3 rows and columns:\n');
fprintf('        Sector1   Sector2   Sector3\n');
for i = 1:3
    fprintf('Sector%d: ', i);
    for j = 1:3
        fprintf('%8.3f  ', table_matrix(i,j));
    end
    fprintf('\n');
end

fprintf('\nColumn sums of TABLE (should be 1):\n');
table_col_sums = sum(table_matrix,1);
for j = 1:12
    fprintf('  Column %2d: %.6f\n', j, table_col_sums(j));
end
