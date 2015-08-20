function ejemplo ();
%--------------------------------------------------------------------------
% ... Este procedimiento carga el archivo de datos *** data_numbers.mat ***
% a la memoria principal. Tambien lo reorganiza en X, una matriz de 
% 400 x 5000. 
% 
% En X, las columnas estan ordenadas en bloques de 500. El primer bloque
% contiene la muestra correspondiente al '0', el segundo bloque contiene la
% muestra del '1', etc. 
%
% curso:      Modelos matematicos y numericos
% instructor: JL Morales
%             Departamento de Matematicas
%             I T A M
%             agosto 2015
%--------------------------------------------------------------------------
%   
load('data_numbers.mat');
X = X';
%
% ... El siguiente segmento de programa elige aleatoriamente una columna
% del X y la despliega como imagen de 20 x 20 en tonos de gris.
%
for i=1:10
    r = randi(5000,1);
    x = X(:,r);
    z = reshape(x,20,20);
    imagesc(z);
    colormap(gray);
    pause
end