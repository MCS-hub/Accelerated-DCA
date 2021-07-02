function x = proximal(y,mu,r)
% compute argmin{ (1/2)||beta - y||^2 + mu*||beta||_1 : ||beta||^2 <= r}

x = soft_thres_opt(y,mu);

if norm(x) >= sqrt(r)
    alpha_plus_one = norm(x)/sqrt(r);
    x = x/alpha_plus_one;
end

end

