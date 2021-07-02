function [obj_list, time, x] = DCA(A,b,sigma,x0,radius,norm_name,tol,omega)

obj_list = [];
obj = objective(A,b,x0,omega);
obj_list = [obj_list,obj];

time = [0];
acc_time = 0;

tic 
I = eye(size(A));
x = x0;
prex = x;

while(1)
    y = (sigma*I-A)*x + omega*sign(x);
    z = -(b-y)/sigma;
    x = project(z,radius,norm_name);

    if norm(prex-x)<tol
        break
    end
    prex = x;
    
    acc_time = acc_time + toc;
    time = [time,acc_time];
    
    obj = objective(A,b,x,omega);
    obj_list = [obj_list,obj];
    
    tic
end
end

