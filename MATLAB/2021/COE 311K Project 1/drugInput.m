function p = drugInput(delta, params, t)

p = 0;

for i = 1:params.Ntreat
    val = delta(i) * (1/(params.sigma*(2*pi)^.5)) * exp(-(abs(params.tau(i)-t)^2)/(2*params.sigma^2));
    p = p + val;
end

end
