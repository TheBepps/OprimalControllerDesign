function k = getMaxK(gamma, k_opt, gamma_opt)

if sort(gamma_opt, 'descend') ~= gamma_opt
    disp('WARNING: the elements of array gamma_opt is not in descending order, result may be imprecise (getMaxK)');
    disp(newline);
end

i=1;
while gamma_opt(i) >= gamma
    i = i+1;
end

if i==1
    k = 0;
else
    k = k_opt(i-1) + (k_opt(i)-k_opt(i-1))*(gamma-gamma_opt(i-1))/(gamma_opt(i)-gamma_opt(i-1));
end