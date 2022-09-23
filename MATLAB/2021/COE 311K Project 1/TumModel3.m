function [g, f, t] = TumModel3(x, params)

deltaT = .01*params.T/(2^params.k); t = linspace(0,params.T,params.T/deltaT+1);

f = zeros(1,length(t)); g = zeros(1,length(t)); f(1) = params.f0; g(1) = params.g0;

for i = 1:length(t) - 1
    
    p = drugInput3(x, params, t(i));
    
    df = -params.Ld*f(i) + p;
    f(i+1) = f(i) + df*deltaT;
    
    dg = params.Lp*g(i)*(1-g(i)) - params.La*g(i) - params.Lk*g(i)*f(i);
    g(i+1) = g(i) + dg*deltaT;

end

end
