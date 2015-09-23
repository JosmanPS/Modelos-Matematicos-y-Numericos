function  ipm_mehrotra_dual();
%--------------------------------------------------------------------------
%... Una prueba del código de puntos interiores extendido a la resolución
% de problemas de programación cuadrática, modificado del original del 
% instructor del curso: 'JL Morales 2013, ITAM'
%
%                                           Omar Pardo 130013
%                                           Modelos Matemáticos
%                                                       2015
%                                                       ITAM
%--------------------------------------------------------------------------

TOL = 1.0e-8;   % ... tolerancia para la brecha de dualidad promedio

% Cargamos los datos
[T,test,ntrain,ntest] = wdbcData('wdbc.data',30,0.0,1);

% Obtenemos las dimensiones y quitamos el identificador
[n_row, n_col] = size(T);
n_atr = 30;
X = T(1:n_row, 2:n_atr+1);

% Escalamos la matriz y la transponemos para dejarla en términos del modelo
[ X ] = scale(X);
X = X';

% Reetiquetamos los valores que tenían 0
ind = T(:,1) ~= 1;
T(ind,1) = -1;

% Definimos las matrices de acuerdo al modelo
b = T(:,1);
Y = diag(b);

A = X*Y; 

gamma = 1000;

fprintf(' Number of samples     ..........  %3i  \n', n_row);
fprintf(' Number of atributes   ..........  %3i  \n', n_atr);

% Le aplicamos puntos interiores con sigma dinámico
[ x, lambda, s ] = ipm_method ( A, b, gamma, TOL );

lambda;

function [ x, y, lambda ] = ipm_method ( A, b, gamma, TOL );
%
%
%--------------------------------------------------------------------------

% Definimos las dimensiones
[n, m] = size(A);

% Definimos el punto inicial para 'x', 's', 'y' y 'lambda'
e = ones(m,1);
x = (gamma/2)*e; 
s = x;
y = A*x;
lambda = 1;

% Calculamos matrices auxiliares
X = diag(x); X = sparse(X);
X_1 = diag(1./x); X_1 = sparse(X_1);
S = diag(s); S = sparse(S);
S_1 = diag(1./s); S_1 = sparse(S_1);

% Definimos z y w iniciales
z = x;
w = x;

% Calculamos matrices auxiliares
Z = diag(z); Z = sparse(Z);
W = diag(w); W = sparse(W);

% Calculamos la brecha inicial
mu = (x'*z+s'*w)/(2*m);

% Iniciamos el vector de F's y definimos la tau
F = zeros(n+m+1);
tau = 0.9995d0;

% Definimos las condiciones de F, sin tomar en cuenta mu

F1 = -e-z-lambda*b+A'*y+w;
F2 = b'*x;
F3 = A*x-y;
F4 = x-gamma*e+s;
F5 = X*z;
F6 = S*w;
F1_bis = F1+X_1*F5-S_1*F6+S_1*W*F4;

% Calculamos las normas de F
F1_n = norm(F1);
F2_n = norm(F2);
F3_n = norm(F3);
F4_n = norm(F4);
F5_n = norm(F5);
F6_n = norm(F6);

% Definimos la función objetivo y la brecha inicial
OBJ =  (0.5)*y'*y-e'*x;
d_gap = mu;

% Calculamos el vector inicial F
F   = -[ F1_bis; F3; -F2 ]; F_n = norm(F); iter = 0;

fprintf('\n');
fprintf('iter   d_gap         OBJ         ||F1||   ||F2||  ||F3||   ||F4||  ||F5||   ||F6||    alpha   sigma\n');
fprintf('------------------------------------\n');
fprintf('%3i   %8.2e  %14.8e %8.2e %8.2e %8.2e %8.2e %8.2e %8.2e \n', ...
        iter, d_gap, OBJ, F1_n, F2_n, F3_n, F4_n, F5_n, F6_n ); 

% Establecemos como condición de paro que la brecha del dual y el primal
% sea menor que la tolerancia y limitamos el método a 30 iteraciones
while d_gap >  TOL & iter < 30 
    
    % Iniciamos el conteo de iteraciones
    iter = iter + 1;    
    
    % Obtenemos la Jacobiana
    KKT = [   X_1*Z + S_1*W         A'         -b     ;
                   A            -eye(n)     zeros(n,1);  
                  -b'           zeros(1,n)     0     ];              
            
    KKT = sparse(KKT);
    
    % Resolvemos el sistema
    [ L, D, P, S ] = ldl(KKT);
   
    dt = backsolve( L, D, P, S, F );
    
    dx =   dt(1:m);
    dy =   dt(m+1:m+n);
    dlambda = dt(m+n+1);
    ds = -F4-dx;
    dz = -X_1*(F5+Z*dx);
    dw = -S_1*(F6+W*ds);
        
    alpha_x = step_d ( x, dx, 1 );
    alpha_s = step_d ( s, ds, 1 );
    alpha_z = step_d ( z, dz, 1 );
    alpha_w = step_d ( w, dw, 1 );
  
    % ... calculamos el parámetro de centrado sigma
    %
    
    mu_aff = ((x + alpha_x*dx)'*(z + alpha_z*dz)+...
              (s + alpha_s*ds)'*(w + alpha_w*dw))/ (2*m);
    sigma  = (mu_aff/mu)^3;
    
    %
    
    F5 = F5+dx.*dz-sigma*mu*e;
    F6 = F6+ds.*dw-sigma*mu*e;
    F1_bis = F1+X_1*F5-S_1*F6+S_1*W*F4;
    
    F   = -[ F1_bis; F3; -F2 ];  
    
     
    % Calculamos el paso corrector
    
    dt = backsolve( L, D, P, S, F );
    
    dx =   dt(1:m);
    dy =   dt(m+1:m+n);
    dlambda = dt(m+n+1);
    
    ds = -F4-dx;
    dz = -X_1*(F5+Z*dx);
    dw = -S_1*(F6+W*ds);
    
    alpha_x = step_d ( x, dx, tau );
    alpha_s = step_d ( s, ds, tau );
    alpha_z = step_d ( z, dz, tau );
    alpha_w = step_d ( w, dw, tau );
    alpha = min(alpha_x, alpha_s);
       %
    % ... movemos las variables según el porcentaje calculado
    %
    x = x + alpha_x*dx;  
    y = y + alpha*dy;
    lambda = lambda + alpha*dlambda;
    s = s + alpha_s*ds;
    z = z + alpha_z*dz;
    w = w + alpha_w*dw;
    
    %
    % ... recalculamos los valores
    %
    X = diag(x); X = sparse(X);
    X_1 = diag(e./x); X_1 = sparse(X_1);
    S = diag(s); S = sparse(S);
    S_1 = diag(e./s); S_1 = sparse(S_1);
    Z = diag(z); Z = sparse(Z);
    W = diag(w); W = sparse(W);
    
    mu = (x'*z+s'*w)/(2*m);       
    d_gap  = mu;
    
    % Redefinimos las condiciones de F
    F1 = -e-z-lambda*b+A'*y+w;
    F2 = b'*x;
    F3 = A*x-y;
    F4 = x-gamma*e+s;
    F5 = X*z;
    F6 = S*w;
    F1_bis = F1+X_1*F5-S_1*F6+S_1*W*F4;

    F   = -[ F1_bis; F3; -F2 ];  F_norm = norm(F);
    
    % Calculamos las normas de F
    F1_n = norm(F1);
    F2_n = norm(F2);
    F3_n = norm(F3);
    F4_n = norm(F4);
    F5_n = norm(F5);
    F6_n = norm(F6);

    % Definimos la función objetivo 
    OBJ =  (0.5)*y'*y-e'*x;
    
    % ... imprimimos los resultados de la iteración anterior
    %
    fprintf('%3i   %8.2e  %14.8e %8.2e %8.2e %8.2e %8.2e %8.2e %8.2e %8.2e %8.2e  \n', ...
        iter, d_gap, OBJ, F1_n, F2_n, F3_n, F4_n, F5_n, F6_n, alpha, sigma );   

end


%
%--------------------------------------------------------------------------
%
function alpha = step_d ( x, dx, tau );
one = 1.0d0; zero = 0.0d0;

ind = find( dx<zero );

if isempty(ind) == 1
    alpha = tau;
else
    alpha = min(-x(ind)./dx(ind));
    alpha = tau*min (one, alpha);
end
%
%--------------------------------------------------------------------------
% ... estandarizar una matriz
%
function [ X ] = scale(X);
[ m, n ] = size(X); e = ones(m,1);

for i=1:n
   xm = mean( X(:,i) );
   xs = std( X(:,i) );
   X(:,i) = (X(:,i) - xm*e)/xs;
end


function [ d ] = backsolve ( L, D, P, S, F );

d = P'*(S*F);
d = L\d;
d = D\d;
d = L'\d;
d = S*P*d;


%
function [train,test,ntrain,ntest] = wdbcData(datafile,dataDim,fracTest,reord)
% syntax: [train,test,ntrain,ntest] = wdbcData(datafile,dataDim,fracTest,reord)
% extract data from the database
% here, "datafile" should be a string, eg 'wdbc.data'
%       "dataDim" is a scalar, eg 30
%       "fracTest" is a scalar strictly between 0 and 1 indicating
%       what fraction of data should go in the test set (the
%       remaining data goes in the training set)
%       "reord" indicates whether the data should be reordered
%       before selecting the test/training sets. Value of "0"
%       indicates no reordering, "1" indicates random reordering
%
% on return, ntrain and ntest indicate the number of rows in the
% training and testing sets, respectively. train and test contain the
% data arrays. Elements in first column of these matrices are either 0
% (indicating a benign sample) or 1 (indicating a malignant
% sample). Elements in the remaining columns, starting with column
% 2, contain the features.


% first check input data
if (nargin < 3)
    error('three or four input arguments are required for wdbcData');
end
if (~ischar(datafile))
    error('first argument must be a string');
end
if (isnumeric(dataDim)&max(size(dataDim))==1)
    if (isempty(dataDim))
        error('second argument must be a scalar');
    elseif (dataDim <= 0)
        error('second argument must be a positive integer');
    end
else
    error('second argument must be a positive number');
end
if(isnumeric(fracTest))
%     if(fracTest<=0 | fracTest>=1)
%         error('third argument must be a number strictly between 0 and 1');
%     end
else
    error('third argument must be a number');
end

if nargin==4
    if(~isnumeric(reord))
        error('fourth argument should be numeric, either 0 or 1');
    end
else
    % default is no reordering
    reord = 0;
end
reord = 0;



fp = fopen(datafile,'r');
samples = 0;
train = zeros(0,dataDim+1);
[id,count] = fscanf(fp,'%d,',1);
while (count == 1)
    samples = samples + 1;
    type = fscanf(fp,'%1s',1);
    if type=='B'
        type =0;
    elseif type=='M'
        type =1;
    else
        type=-1;
        fprintf(' invalid type found in data file %s\n', datafile);
    end
    vec = fscanf(fp,',%e',dataDim);
    train = [train; type vec'];
    [id,count] = fscanf(fp,'%d,',1);
end
if (samples < 569)
    error('Not enough samples');
end

% reorder the rows of "train", if requested
if ~(reord==0)
    p = randperm(samples);
    train(p,:) = train;
end

ntest = round(fracTest * samples);
if ntest < 1
    ntest=1;
elseif ntest >= samples
    ntest = samples-1;
end
ntrain = samples - ntest;

% test = train(1:ntest,:);
% train = train(ntest+1:samples,:);
% modified 4/7/03
test  = train(ntrain+1:samples,:);
train = train(1:ntrain,:);

%--------------------------------------------------------------------------
%