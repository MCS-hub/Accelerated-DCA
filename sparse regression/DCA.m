function [time,objective_list,beta] = DCA(X,y,beta0,alpha,lambda,sigma,r,tol)
disp('DCA')
objective_list = [];
beta = beta0;
objective = objective_value(beta,X,y,lambda,alpha);
objective_list = [objective_list,objective];

time = [0];
acc_time = 0;

tic
A = X'*X;
b = X'*y;

pre_beta = beta;
while true
    disp('count')
    v = gradH(beta,A,alpha,lambda,sigma);
    z = (b+v/2)/sigma;
    mu = (lambda*alpha)/(2*sigma);
    beta = proximal(z,mu,r);
    
    acc_time = acc_time + toc;
    time = [time, acc_time];
    
    objective = objective_value(beta,X,y,lambda,alpha);
    objective_list = [objective_list,objective];
    
    tic
    
    if norm(beta-pre_beta)<tol
        break
    end
    
    pre_beta = beta;
end
end

