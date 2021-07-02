function v = gradH(beta,A,alpha,lambda,sigma)

v = 2*sigma*beta - 2*A*beta + lambda*alpha*sign(beta).*(alpha*abs(beta)>1);

end

