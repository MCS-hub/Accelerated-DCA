function objective = objective_value(beta,X,y,lambda,alpha)

objective = norm(y-X*beta)^2 + lambda*sum(min(1,alpha*abs(beta)));

end

