function x_proj = project(x,radius,norm_name)

if norm_name == "l2"
    if norm(x) <= radius
        x_proj = x;
    else
        x_proj = radius*x/norm(x);
    end
elseif norm_name == "linf"
    x(x>radius) = radius;
    x(x<-radius) = -radius;
    x_proj = x;
end
end

