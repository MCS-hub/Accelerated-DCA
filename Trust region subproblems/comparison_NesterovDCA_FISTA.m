%
N = 1000;
radius_list = [200];
scale_sigma = [1,5,10,25];
n_runs = 20;

result_path = '.\new_results\';
for i_count = 1:length(radius_list)
    radius = radius_list(i_count);
    for j_count = 1:length(scale_sigma)
        scale = scale_sigma(j_count);
        for k_count = 1:n_runs
            data = load(strcat(result_path,string(radius),'_',string(scale),'_',string(k_count),'.mat'));
            A = data.A;
            b = data.b;
            x0 = data.x0;
            
            norm_name = "linf";
            tol = 1e-6;

            A_max_eig = eigs(A,1,'largestreal');
            sigma = scale*max(1,A_max_eig)+1;
            
            disp('NesDCA')
            q = 5;
            [obj_list_NesDCA, time_NesDCA, x_NesDCA] = ADCA(A,b,sigma,x0,radius,norm_name,tol,q);

            disp('FISTA')
            [obj_list_fista,time_fista,x_fista] = FISTA(A,b,sigma,x0,radius,norm_name,tol);
            
            
            obj_list_DCA = data.obj_list_DCA;
            obj_list_aDCA = data.obj_list_aDCA;
            time_DCA = data.time_DCA;
            time_aDCA = data.time_aDCA;
            
            figure
            min_value = min([obj_list_DCA,obj_list_aDCA,obj_list_NesDCA,obj_list_fista]);
            plot(time_DCA, obj_list_DCA-min_value,'LineWidth',2)
            hold on
            plot(time_aDCA, obj_list_aDCA-min_value,'--','LineWidth',2)
            hold on
            plot(time_NesDCA,obj_list_NesDCA-min_value,'LineWidth',3)
            hold on
            plot(time_fista,obj_list_fista-min_value,'-.','LineWidth',2)
            
            xlabel('time (s)')
            ylabel('f(x_k) - f^*')
            legend('DCA','accelerated DCA','Nesterov DCA','FISTA')
            set(gca,'YScale','log')
            
            savefig(strcat(result_path,string(radius),'_',string(scale),'_',string(k_count),'.fig'));
            
            save(strcat(result_path,string(radius),'_',string(scale),'_',string(k_count),'.mat'),'A','b','x0','obj_list_DCA',...
                'time_DCA','x_DCA','obj_list_aDCA','time_aDCA','x_aDCA','time_NesDCA','x_NesDCA','obj_list_NesDCA',...
                'time_fista','x_fista','obj_list_fista');
            close all
        end
    end
end
