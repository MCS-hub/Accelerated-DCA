function [time,objective_list,beta] = accelerated_DCA(X,y,beta0,alpha,lambda,sigma,r,tol,ls_armijo_bar,gamma,rho,bst_dis)
disp('ADCA+')
[N,m] = size(X);

objective_list = [];
beta = beta0;
objective = objective_value(beta,X,y,lambda,alpha);
objective_list = [objective_list,objective];

time = [0];
acc_time = 0;

tic
A = X'*X;
b = X'*y;

while(1)
    v = gradH(beta,A,alpha,lambda,sigma);
    z = (b+v/2)/sigma;
    mu = (lambda*alpha)/(2*sigma);
    beta_tmp = proximal(z,mu,r);
    
    d = beta_tmp - beta;  % boosting direction
    ns_d = norm(d)^2;
    
    disp('norm beta square')
    disp(norm(beta_tmp)^2)
    
    if ns_d < tol^2
        break
    end
    
    if ns_d < bst_dis
        disp('boosting...')
        sg_h_1d = 2*sigma*d'*beta_tmp - 2*d'*A*beta_tmp +...
            lambda*alpha*sum(d.*sign(beta_tmp).*(1<=alpha*abs(beta_tmp)));
        
        tmp1 = d'*beta_tmp;
        tmp2 = sqrt(ns_d*(r-norm(beta_tmp)^2)+tmp1^2);
        s2 =(tmp2-tmp1)/ns_d;
        s1 = (-tmp2-tmp1)/ns_d;
        
        tmp3 = 2*tmp1*sigma - 2*d'*b - sg_h_1d;
        
        % solve the 1D convex problem using CPlex
        H = zeros(m+1);
        H(1,1) = 2*sigma*ns_d;
        f = [tmp3;lambda*alpha*ones(m,1)];
        Aineq = [d,-eye(m);-d,-eye(m)];
        bineq = [-beta_tmp;beta_tmp];
        lb = [s1;-inf*ones(m,1)];
        ub = [s2; inf*ones(m,1)];
        sol = cplexqp(H,f,Aineq,bineq,[],[],lb,ub);
        s = sol(1);
        
        disp('s')
        disp(s)
        if s>0
            % Armijo-type line-search
            ls_armijo = ls_armijo_bar;
            
            %objective(A,b,v+lambda*d) > objective(A,b,v)-(rho/2)*ns_d*lambda^2 ...
            obj_beta_tmp = objective_value(beta_tmp,X,y,lambda,alpha);
            
            while objective_value(beta_tmp+ls_armijo*d,X,y,lambda,alpha) > obj_beta_tmp...
                    -(rho/2)*ns_d*ls_armijo^2 || ls_armijo < s1 || ls_armijo > s2
                ls_armijo = gamma*ls_armijo;
            end
            disp('ls_armijo')
            disp(ls_armijo)
        else
            ls_armijo = s;
        end
        
        beta = beta_tmp + ls_armijo*d;
        
    else
        beta = beta_tmp;
    end
    
    acc_time = acc_time + toc;
    time = [time, acc_time];
    
    objective = objective_value(beta,X,y,lambda,alpha);
    objective_list = [objective_list,objective];
    
    tic    
    
end
end

