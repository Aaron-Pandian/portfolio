function error = drugError(delta, params)

finalG = zeros(1,11);
finalG(1)=params.g0;

for k = 0:10
    params.k = k;
    [g,f,t] = TumModel(delta, params);
    finalG(k+1) = g(end);
end


error = zeros(1,10);

for k0 = 1:10
    error(k0) = 100*(finalG(k0+1) - finalG(k0))/finalG(k0);
end