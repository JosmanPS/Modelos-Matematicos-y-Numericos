classdef digits
        
% DOC
    
    properties (SetObservable = true)
        k = 2       % Número de valores singulares usados
        tasa        % Tasa de acierto
        x           % Vector de datos x
        y_real      % Vector respuesta
        y_pred      % Predicción
        flag        % Acierto
    end

    properties (SetAccess = private, Hidden = true)
        U_one
        U_two
        U_three
        U_four
        U_five
        U_six
        U_seven
        U_eight
        U_nine
        U_zero
    end
    
  
    methods
        function obj = init(obj)
            u = load('U_Matrix.mat');
            obj.U_one = u.U_one;
            obj.U_two = u.U_two;
            obj.U_three = u.U_three;
            obj.U_four = u.U_four;
            obj.U_five = u.U_five;
            obj.U_six = u.U_six;
            obj.U_seven = u.U_seven;
            obj.U_eight = u.U_eight;
            obj.U_nine = u.U_nine;
            obj.U_zero = u.U_zero;
        end

        function obj = init_ml(obj)
            u = load('test_matrix.mat');
            obj.U_one = u.U_one;
            obj.U_two = u.U_two;
            obj.U_three = u.U_three;
            obj.U_four = u.U_four;
            obj.U_five = u.U_five;
            obj.U_six = u.U_six;
            obj.U_seven = u.U_seven;
            obj.U_eight = u.U_eight;
            obj.U_nine = u.U_nine;
            obj.U_zero = u.U_zero;
        end

        function obj = random(obj)
            load('data_numbers.mat');
            index = randi([1 5000]);
            x = X(index, :);
            obj.x = x';
            obj.y_real = y(index);
            clear X
            clear y
        end

        function obj = random_ml(obj)
            load('test.mat');
            index = randi([1 1000]);
            x = X_test(index, :);
            obj.x = x';
            obj.y_real = y_test(index);
            clear X
            clear y
        end

        function obj = predict(obj)

            x = obj.x;
            k = obj.k;

            r = zeros(10, 1);

            r(1) = norm( x - obj.U_one(:, 1:k) * obj.U_one(:, 1:k)' * x );
            r(2) = norm( x - obj.U_two(:, 1:k) * obj.U_two(:, 1:k)' * x );
            r(3) = norm( x - obj.U_three(:, 1:k) * obj.U_three(:, 1:k)' * x );
            r(4) = norm( x - obj.U_four(:, 1:k) * obj.U_four(:, 1:k)' * x );
            r(5) = norm( x - obj.U_five(:, 1:k) * obj.U_five(:, 1:k)' * x );
            r(6) = norm( x - obj.U_six(:, 1:k) * obj.U_six(:, 1:k)' * x );
            r(7) = norm( x - obj.U_seven(:, 1:k) * obj.U_seven(:, 1:k)' * x );
            r(8) = norm( x - obj.U_eight(:, 1:k) * obj.U_eight(:, 1:k)' * x );
            r(9) = norm( x - obj.U_nine(:, 1:k) * obj.U_nine(:, 1:k)' * x );
            r(10) = norm( x - obj.U_zero(:, 1:k) * obj.U_zero(:, 1:k)' * x );

            [ ~, obj.y_pred ] = min(r);
            obj.flag = obj.y_real == obj.y_pred;

        end

        function tasa = run(obj, n)
            results = zeros(n, 1);
            
            for i = 1:n
                o = random(obj);
                o = predict(o);
                results(i) = o.flag;
            end
            
            tasa = mean(results);

        end

        function tasa = run_ml(obj, n)
            results = zeros(n, 1);
            
            for i = 1:n
                o = random_ml(obj);
                o = predict(o);
                results(i) = o.flag;
            end
            
            tasa = mean(results);

        end
        
    end
    
end
