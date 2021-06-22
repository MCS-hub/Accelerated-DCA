%
N = 1000;
radius_list = [200,500,1000];
scale_sigma = [1,5,10,25];
n_runs = 20;
result_path = '.\results\';
for i_count = 1:length(radius_list)
    radius = radius_list(i_count);
    for j_count = 1:length(scale_sigma)
        scale = scale_sigma(j_count);
        for k_count = 1:n_runs
            A = triu(randn(N,N));
            B = A';
            A = A+B;
            b = randn(N,1);

            radius = 200;
            x0 = randn(N,1);
            norm_name = "linf";
            tol = 1e-8;

            A_max_eig = eigs(A,1,'largestreal');
            sigma = scale*max(1,A_max_eig)+1;

            disp('DCA')
            [obj_list_DCA, time_DCA, x_DCA] = DCA(A,b,sigma,x0,radius,norm_name,tol);

            disp('accelerated DCA')
            %rho = sigma - A_max_eig;
            rho = 1;
            lambda0 = 0.01;
            gama = 0.1; % 0.5?
            xi = 5.;
            [obj_list_aDCA, time_aDCA, x_aDCA] = acceleratedDCA(A,b,sigma,x0,radius,norm_name,tol,rho,lambda0,gama,xi);
            
            figure
            min_value = min([obj_list_DCA,obj_list_aDCA]);
            plot(time_DCA, obj_list_DCA-min_value,'LineWidth',2)
            hold on
            plot(time_aDCA, obj_list_aDCA-min_value,'--','LineWidth',2)
            xlabel('time (s)')
            ylabel('f(x_k) - f^*')
            legend('DCA','accelerated DCA')
            set(gca,'YScale','log')
            
            savefig(strcat(result_path,string(radius),'_',string(scale),'_',string(k_count),'.fig'));
            
            save(strcat(result_path,string(radius),'_',string(scale),'_',string(k_count),'.mat'),'A','b','x0','obj_list_DCA',...
                'time_DCA','x_DCA','obj_list_aDCA','time_aDCA','x_aDCA');
            close all
        end
    end
end
