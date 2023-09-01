%% DPI
x =[-2 0 2];
y =[4 2 8];
f = polyfit(x,y,2);
disp(f);
t= @(x)(x.^2+x+2);
x_values = -4:0.25:4;
plot(x_values,t(x_values))

%% NPI
x2 = NIP(x,y,-4);
y2 = NIP(x,y,4);
xnew = [-4 -2 0 2 4];
ynew = [q 4 2 8 y2];
fnew = polyfit(xnew,ynew,2);
disp(fnew);
plot(x_values,t(x_values))

%% Lagrange  
x3 = NIP(x,y,-4);
y3 = Lagrange(x,y,4);
xnew2 = [-4 -2 0 2 4];
ynew2 = [w 4 2 8 y3];
fnew2 = polyfit(xnew2,ynew2,2);
disp(fnew2);




