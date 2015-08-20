
% --------------------------------------------------
%
%                 OOP run
% 
% Este script clasifica 100 dígitos con el clasificador
% construido con 1, 2, ..., 20 valores singulares para
% medir su tasa de acierto.
%
% José Manuel Proudinat Silva
% 000130056
%
% --------------------------------------------------

% Creamos el clasificador
clasificador = digits();
clasificador = init(clasificador);

% Clasificamos y guardamos los resultados
tasas_train = zeros(10, 1);
tasas_test = zeros(10, 1);

for k = 1:20
    clasificador.k = k;
    tasas_train(k) = run(clasificador, 100);
end

clasificador = init_ml(clasificador);

for k = 1:20
    clasificador.k = k;
    tasas_test(k) = run_ml(clasificador, 100);
end


plot(1:20, tasas_train, 1:20, tasas_test, 'k', 'LineWidth', 1);
title('Tasa de acierto de entrenamiento (negra) y de prueba (azul)');
xlabel('Número de valores singulares');
ylabel('Tasa de acierto');

