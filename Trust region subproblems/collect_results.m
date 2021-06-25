result_path = '.\results\';
radius = 200;
scale_sigma = [1,5,10,25];
n_runs = 20;

final_time = [];
final_obj = [];
for j_count = 1:length(scale_sigma)
    scale = scale_sigma(j_count);
    for k_count = 1:n_runs
        data = load(strcat(result_path,string(radius),'_',string(scale),'_',string(k_count),'.mat'));
        obj_list_DCA = data.obj_list_DCA;
        time_DCA = data.time_DCA;
        obj_list_aDCA = data.obj_list_aDCA;
        time_aDCA = data.time_aDCA;
        
        obj_DCA_end = obj_list_DCA(end);
        time_DCA_end = time_DCA(end);
        obj_aDCA_end = obj_list_aDCA(end);
        time_aDCA_end = time_aDCA(end);
        
        final_time = [final_time; [time_DCA_end,time_aDCA_end]];
        final_obj = [final_obj; [obj_DCA_end,obj_aDCA_end]];
    end
end

%%
scale = 25;
bar(final_time(61:80,:))
xlabel('dataset')
ylabel('time (s)')
legend('DCA','accelerated DCA')
savefig(strcat('time_scaleSigma_',string(scale),'.fig'));

bar(final_time(21:40,:))
bar(final_time(41:60,:))
bar(final_time(61:80,:))
%%
x_axis = 1:80;
plot(x_axis,final_obj(:,1),'b','LineWidth',2.5);
hold on
plot(x_axis,final_obj(:,2),'r--','LineWidth',2);
legend('DCA','accelerated DCA');
xlabel('dataset');
ylabel('final objective value');