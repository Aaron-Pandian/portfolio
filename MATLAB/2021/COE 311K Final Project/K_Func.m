function k = K_Func(x, k1, k2, L1, L12)
if x <= L1
    k = k1;
elseif x >= L1 + L12
    k = k2;
else
    k = k1 + (k2 - k1)*(x - L1)/(L12);
end
end
