% Quick diagnostic: what are the date cells in pib_sectorial_bc.xlsx?
dp  = 'C:\Users\adiaz\OneDrive - Banco Central de Chile\Proyectos\Networks\data_github';
f   = fullfile(dp, 'pib_sectorial_bc.xlsx');

% Read just first 10 rows of col A (row 3 = header, rows 4-12 = first data)
raw = readcell(f, 'Sheet', 'Cuadro', 'Range', 'A3:A12');
fprintf('Total rows read: %d\n', size(raw,1));
for k = 1:size(raw,1)
    d = raw{k,1};
    cl = class(d);
    if isdatetime(d)
        fprintf('Row %d [%s]: datetime = %s  (year=%d, month=%d)\n', k, cl, char(d), year(d), month(d));
    elseif isnumeric(d) && ~isnan(d)
        fprintf('Row %d [%s]: numeric = %g\n', k, cl, d);
    elseif ischar(d) || isstring(d)
        fprintf('Row %d [%s]: string = "%s"\n', k, cl, char(d));
    elseif ismissing(d)
        fprintf('Row %d [%s]: MISSING\n', k, cl);
    else
        fprintf('Row %d [%s]: other\n', k, cl);
    end
end

% Also check what the full row 4 (first data row) looks like for cols A-C
raw2 = readcell(f, 'Sheet', 'Cuadro', 'Range', 'A4:C6');
fprintf('\nFirst 3 data rows, cols A-C:\n');
for r = 1:size(raw2,1)
    for c = 1:size(raw2,2)
        d = raw2{r,c};
        fprintf('  [%d,%d] class=%s', r, c, class(d));
        if isnumeric(d) && ~isnan(d)
            fprintf(' val=%g', d);
        elseif ischar(d)||isstring(d)
            fprintf(' val="%s"', char(d));
        elseif isdatetime(d)
            fprintf(' val=%s', char(d));
        end
        fprintf('\n');
    end
end

% Also list sheet names
[~, snames] = xlsfinfo(f);
fprintf('\nSheet names in file:\n');
for k = 1:numel(snames)
    fprintf('  %s\n', snames{k});
end
