[g1, f1, g2, f2, t]=solveOptimization();
[g3, f3]=solveOptimization3();

subplot(2,1,1)
plot(t, g1)
hold on
plot(t, g2)
hold on
plot (t, g3)
legend('initial g','optimal g', 'new optimal g')

subplot(2,1,2)
plot(t, f1)
hold on 
plot(t, f2)
hold on
plot(t, f3)
hold on
legend('initial f', 'optimal f', 'new optimal f')
