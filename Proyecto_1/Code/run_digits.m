function [ tasa ] = run_digits( n, k )

% --------------------------------------------------
%
%                     RUN DIGITS
% 
% Esta función usa el predictor de dígitos un número
% definido de veces y devuelve la tasa de acierto
% obtenida.
%
% INPUT:
%     - n     : El número de veces que se evaluará
%               el predictor.
%     - k     : El número de valores singulares que
%               se utilizarán para la predicción.
%
% OUTPUT:
%     - tasa : La tasa de acierto del predictor.
%
% José Manuel Proudinat Silva
% 000130056
%
% --------------------------------------------------

    % Creamos el vector que almacenará los resultados
    % del predictor
    results = zeros(n, 1);
    
    for i = 1:n
        [ ~, r ] = digit_class( k );
        results(i) = r;
    end

    tasa = mean(results);
    
end
        