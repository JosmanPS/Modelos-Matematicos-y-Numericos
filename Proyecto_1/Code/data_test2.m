clear all

load('data_numbers.mat')

for i = 1:10
    X_test((1 + 100*(i-1)):(100*i), :) = X((401 + 500*(i-1)):(500*i), :);
end

base = ones(1, 100);
y_test = [10*base, base, 2*base, 3*base, 4*base, 5*base, 6*base, 7*base, 8*base, ...
         9*base];

clear X
clear y

save('test.mat')