function x = proximal(y,mu,r)
% compute argmin{ (1/2)||x - y||^2 + mu*||x||_1 : ||x||^2 <= r}

x = soft_thres_opt(y,mu);
norm_x = norm(x);
if norm_x >= sqrt(r)
    alpha_plus_one = norm_x/sqrt(r);
    x = x/alpha_plus_one;
end

end

