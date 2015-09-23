
[T,test,ntrain,ntest] = wdbc('wdbc.data',30,0.0,1);

[n_row, n_col] = size(T);
n_atr = 30;
X = T(1:n_row,2:n_atr+1);
X = scale(X);
% Reetiquetamos los valores que tenían 0
ind = T(:,1) ~= 1;
T(ind,1) = -1;
y = T(:, 1);

mu = 1.0e3;     % ... penalizacion para el termino lineal
TOL = 1.0e-8;   % ... tolerancia para la brecha de dualidad promedio

for i=1:n_row
   if T(i,1) == 1
       m = m + 1;
       M(m,:) = T(i,2:31);
   else
       k = k + 1;
       B(k,:) = T(i,2:31);
   end
end

mm = m + k;
nn = 2*n_atr + 2 + m + k + m + k;

fprintf(' Number of samples     ..........  %3i  \n', n_row);
fprintf(' Number of atributes   ..........  %3i  \n', n_atr);
fprintf(' Number of malignant samples ....  %3i  \n', m);
fprintf(' Number of benign samples .......  %3i  \n', k);
fprintf(' Number of variables      .......  %3i  \n', nn);
fprintf(' Number of constraints    .......  %3i  \n', mm);

[ x, y, s ] = svm_dual_mehrotra(X, y, mu, TOL);