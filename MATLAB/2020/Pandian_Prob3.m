x = randn(10^5,1);
figure;
histogram(x);
old_size = length(x);
A = x > -1;
B = 1 < x;
new_x = A | B;
new_x(new_x==0) = [];   
new_size = length(new_x);
percentage_removed = 100 * ((old_size - new_size) / (new_size));
disp(percentage_removed)