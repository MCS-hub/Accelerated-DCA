N = 1000;
radius = 200;
scale = 10;

A = triu(randn(N,N));
B = A';
A = A+B;
b = randn(N,1);

x0 = randn(N,1);
norm_name = "linf";
tol = 1e-6;
omega = 100;

A_max_eig = eigs(A,1,'largestreal');
sigma = scale*max(1,A_max_eig)+1;

[obj_list_DCA, time_DCA, x_DCA] = DCA(A,b,sigma,x0,radius,norm_name,tol,omega);

bst_dis = 0.1;
rho = 1;
lambda0 = 100;
gama = 0.5;
[obj_list_aDCA,time_aDCA,x_aDCA] = acceleratedDCA_v2(A,b,sigma,x0,radius,norm_name,tol,rho,lambda0,gama,bst_dis,omega);

figure
min_value = min([obj_list_DCA,obj_list_aDCA]);
plot(time_DCA, obj_list_DCA-min_value,'LineWidth',2)
hold on
plot(time_aDCA, obj_list_aDCA-min_value,'--','LineWidth',2)
xlabel('time (s)')
ylabel('f(x_k) - f^*')
legend('DCA','accelerated DCA')
set(gca,'YScale','log')