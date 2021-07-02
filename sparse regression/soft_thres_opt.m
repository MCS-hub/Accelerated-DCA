function x = soft_thres_opt(y,mu)

x = sign(y).*max(abs(y)-mu,0);

end

