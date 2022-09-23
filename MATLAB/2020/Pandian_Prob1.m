%% Linear Equations
A = [4,3;-3,2];
disp('A: ');
disp(A)

B = [5;2];
disp('B: ');
disp(B)

fprintf('This is like %dx + %dy = %d \n', A(1,1), A(1,2), B(1));
fprintf('And %dx + %dy = %d \n', A(2,1), A(2,2), B(2));
disp('If cx + dy = b, dy = b - cx, and y = b/d - x(c/d)');
cCol = A(:,1);
dCol = A(:,2);

x = linspace(-3,3,5);
y = (B./dCol - (cCol./dCol).*x);

figure(1);
plot(x,y);
hold on;
solution = A\B;
plot(solution(1),solution(2),'gs','Markersize',10);


