function [SSI] = measure_SSI(A, B)

numrow = size(A, 2);
numcol = size(B, 2);

%%
covAB = zeros(numrow, numcol);
for i = 1:numrow
    for j = 1:numcol
         sigma = corrcoef(A(:, i), B(:, j));
         covAB(i, j) = abs(sigma(1, 2));
         if isnan(covAB(i,j))
             covAB(i, j) = 0;
         end
    end
end

pos_row = find(sum(covAB, 2)>0);
pos_col = find(sum(covAB, 1)>0);

covAB = covAB(pos_row, pos_col);

numrow = numel(pos_row);
numcol = numel(pos_col);

%%
SSI_row = 0;
[max_val_row, max_row] = max(covAB, [], 2);
mean_row = mean(covAB, 2);
SSI_row = SSI_row + sum(max_val_row);
for i = 1:numrow
    ind = setdiff( find(covAB(i, :) > mean_row(i)), max_row(i));
    SSI_row = SSI_row - sum(covAB(i, ind))/(numcol-1 + eps);
end

SSI_row = SSI_row/(2*numrow);

%%
SSI_col = 0;
[max_val_col, max_col] = max(covAB, [], 1);
mean_col = mean(covAB, 1);
SSI_col = SSI_col + sum(max_val_col);
for j = 1:numcol
    ind = setdiff( find(covAB(:, j) > mean_col(j)), max_col(j));
    SSI_col = SSI_col - sum(covAB(ind, j))/(numrow-1 + eps);
end

SSI_col = SSI_col/(2*numcol);

SSI = SSI_col + SSI_row;

if isnan(SSI)
    SSI = 0;
else
    return;
end

end