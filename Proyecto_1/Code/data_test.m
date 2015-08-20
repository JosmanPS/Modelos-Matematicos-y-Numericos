clear all

% Cargamos los datos originales
load('data_numbers.mat')

% Guardamos las matrices de datos
X_test = zeros(1000, 400);

X_zero = X(1:400, :);;
X_one = X(501:900, :);
X_two = X(1001:1400, :);
X_three = X(1501:1900, :);
X_four = X(2001:2400, :);
X_five = X(2501:2900, :);
X_six = X(3001:3400, :);
X_seven = X(3501:3900, :);
X_eigth = X(4001:4400, :);
X_nine = X(4501:4900, :);

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
clear X_one
clear X_two
clear X_three
clear X_four
clear X_five
clear X_six
clear X_seven
clear X_eight
clear X_nine
clear X_zero

% Exportamos los datos
save('test_matrix.mat')