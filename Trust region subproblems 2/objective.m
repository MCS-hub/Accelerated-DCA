function obj = objective(A,b,x,omega)

obj = (1/2)*x'*A*x + b'*x - omega*norm(x,1);

end

