% data_path = '.\standardized_data\';
% data_name = 'eunite2001';
% 
% data = load(strcat(data_path,data_name,'.mat'));
% 
% X = data.X;
% y = data.y;
% [N,d] = size(X);
N = 10000;
m = 1000;
X = randn(N,m);
beta_init = randn(m,1);
noise = 0.05*randn(N,1);
y = X*beta_init + noise;

r = 10000000;
alpha = 1;
lambda = 0.5*N;

A = X'*X;

A_max_eig = eigs(A,1,'largestreal');
sigma = max(1,A_max_eig)+1;
tol = 1e-6;
beta0 = randn(m,1);

[time,objective_list,beta] = DCA(X,y,beta0,alpha,lambda,sigma,r,tol);

ls_armijo_bar = 10;
gamma = 0.5;
rho = 1;
bst_dis = 0.01; %0.001;
[time2,objective_list2,beta2] = accelerated_DCA(X,y,beta0,alpha,lambda,sigma,r,tol,ls_armijo_bar,gamma,rho,bst_dis);

min_value = min([objective_list,objective_list2]);
plot(time,objective_list-min_value)
hold on
plot(time2,objective_list2 - min_value,'--r')
legend('DCA','accelerated DCA')
set(gca,'YScale','log')