function [obj_list,time,x] = FISTA(A,b,sigma,x0,radius,norm_name,tol)

% moi chi cai dat cho L_inf
if norm_name ~= "linf"
    error('invalid norm_name, this norm is not available for a moment...')
end

time = [0];
acc_time = 0;

obj_list = [];
obj = objective(A,b,x0);
obj_list = [obj_list,obj];

tic
I = eye(size(A));
x = x0;

t = 1;

first_iter = true;

while(1)
    t = (1+sqrt(1+4*t^2))/2;
    if first_iter
        z = x0;
        first_iter = false;
    else
        z = x + (pret-1)/t*(x-prex);
    end
    prex = x;
    pret = t;
    
    y = (sigma*I - A)*z - b;
    
    x = project(y/sigma,radius,norm_name);
    
    if norm(prex-x)<tol
        break
    end
    
    acc_time = acc_time + toc;
    time = [time,acc_time];
    
    % compute the objective
    obj = objective(A,b,x);
    obj_list = [obj_list,obj];
    tic
    
end
end