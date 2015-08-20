function [ y , flag ] = digit_class( k )

% --------------------------------------------------
%
%               DIGIT CLASSIFICATOR
% 
% Esta función elige una observación aleatoria de la
% matriz de datos y luego clasifica según el residuo
% que deja la proyección de sus datos en la submatriz
% U_K de la DVS de los datos de entrenamiento de 
% cada dígito.
%
% INPUT:
%     - k     : El número de valores singulares que se
%               utilizarán para la predicción (k <= 20)
%
% OUTPUT:
%     - y     : El dígito que se predijo
%     - flag  : Booleano que representa si se tuvo
%               éxito o error en la clasificación
%
% José Manuel Proudinat Silva
% 000130056
%
% --------------------------------------------------

    % Obtenemos la observación aleatoria
    index = randi([1 5000]);
    load('data_numbers.mat');
    x = X(index, :);
    x = x';
    y_real = y(index);
    clear X
    clear y

    % Cargamos las matrices U
    load('U_Matrix.mat');

    % Verificamos el residuo para la proyección en cada
    % dígito y lo almacenamos
    r = zeros(10, 1);

    r(1) = norm( x - U_one(:, 1:k) * U_one(:, 1:k)' * x );
    r(2) = norm( x - U_two(:, 1:k) * U_two(:, 1:k)' * x );
    r(3) = norm( x - U_three(:, 1:k) * U_three(:, 1:k)' * x );
    r(4) = norm( x - U_four(:, 1:k) * U_four(:, 1:k)' * x );
    r(5) = norm( x - U_five(:, 1:k) * U_five(:, 1:k)' * x );
    r(6) = norm( x - U_six(:, 1:k) * U_six(:, 1:k)' * x );
    r(7) = norm( x - U_seven(:, 1:k) * U_seven(:, 1:k)' * x );
    r(8) = norm( x - U_eight(:, 1:k) * U_eight(:, 1:k)' * x );
    r(9) = norm( x - U_nine(:, 1:k) * U_nine(:, 1:k)' * x );
    r(10) = norm( x - U_zero(:, 1:k) * U_zero(:, 1:k)' * x );

    % Observamos el mínimo
    [~, y] = min(r);
    flag = y == y_real;

    if ( y == 10)
       y = 0;
    end

end


    