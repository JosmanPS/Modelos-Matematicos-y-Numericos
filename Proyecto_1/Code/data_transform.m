
% --------------------------------------------------
%
%                 DATA TRANSFORM
% 
% Este script toma las matrices de datos originales
% y las separa según el dígito que representan para
% hacer más intuitivo el código posterior. Además
% resuelve el problema DVS para cada una y exporta
% las U_i necesarias para la clasificación. (K=20)
%
% José Manuel Proudinat Silva
% 000130056
%
% --------------------------------------------------

clear all

% Cargamos los datos originales
load('data_numbers.mat')

% Guardamos las matrices de datos
X_zero = X(1:500, :);
X_one = X(501:1000, :);
X_two = X(1001:1500, :);
X_three = X(1501:2000, :);
X_four = X(2001:2500, :);
X_five = X(2501:3000, :);
X_six = X(3001:3500, :);
X_seven = X(3501:4000, :);
X_eigth = X(4001:4500, :);
X_nine = X(4501:5000, :);

clear X
clear y

% Calculamos la DVS y guardamos la U_K correspondiente
% Además quitamos de la memoria lo que no necesitamos
[U, ~, ~] = svd(X_zero');
U_zero = U(:, 1:20);
[U, ~, ~] = svd(X_one');
U_one = U(:, 1:20);
[U, ~, ~] = svd(X_two');
U_two = U(:, 1:20);
[U, ~, ~] = svd(X_three');
U_three = U(:, 1:20);
[U, ~, ~] = svd(X_four');
U_four = U(:, 1:20);
[U, ~, ~] = svd(X_five');
U_five = U(:, 1:20);
[U, ~, ~] = svd(X_six');
U_six = U(:, 1:20);
[U, ~, ~] = svd(X_seven');
U_seven = U(:, 1:20);
[U, ~, ~] = svd(X_eigth');
U_eight = U(:, 1:20);
[U, ~, ~] = svd(X_nine');
U_nine = U(:, 1:20);

clear U
clear x

% Exportamos los datos
save('U_Matrix.mat')
