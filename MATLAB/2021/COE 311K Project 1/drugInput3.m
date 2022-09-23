function p = drugInput3(x, params, t)

p = 0;
delta=x(1:4);
tau=x(5:8);


for i = 1:params.Ntreat
    val = delta(i) * (1/(params.sigma*(2*pi)^.5)) * exp(-(abs(tau(i)-t)^2)/(2*params.sigma^2));
    p = p + val;
end
end
