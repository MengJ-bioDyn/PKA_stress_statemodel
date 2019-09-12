%%Code for plot any single trajectory from the data_final_s2

%%
clear all

load data_final_s2;
%Type in the index of cells need to be plotted
selected = input('please type in the cell you want to plot ...[1 2] :   ');

for i=1:length(all_data_2)
    x(i)=length(all_data_2(i).traj_i);
    y(i)=max(all_data_2(i).traj_i);
    
    date = all_data(i).Date;
        date = num2str(date);
        if strcmp(date,'20151112')
            interval=6;
        else
            interval=15;
        end
    
        %Get the cell cycle info
        cycles = all_data(i).cycle(5:end);
        %Remove the zeros
        cycles = cycles(cycles>0);
        cycle_start=all_data(i).cycle(3);
        
        all_data_2(i).cycles=(cycles-cycle_start)*interval;
        
        all_data_2(i).age=all_data(i).age;
end
% figure;
% plot((1:length(all_data_2(1).traj_s)),all_data_2(1).traj_s,'r-','LineWidth',2);
% y1 = get(gca,'ylim');
%         for k_y = all_data_2(1).cycles
%             line([k_y k_y],y1,'Color','k','LineStyle','--','LineWidth',1.5)
%         end


x_max=max(x);
y_max=max(y);

for CCI = selected

figure


% Initiating subplot
% For each figure

    
   
    



plot((1:length(all_data_2(CCI).traj_s)),all_data_2(CCI).traj_s,'r-','LineWidth',2);
y1 = get(gca,'ylim');
y1 = [0 4];
% hold on;
% plot((1:length(all_data_2(CCI).traj_i)),all_data_2(CCI).traj_i,'b.','MarkerSize',5);
ax1 = gca;
set(ax1,'FontSize',10,'FontName','Arial')
xlabel('Time (min)');
%,'XTick',[0.01 0.1 1 10],'YTick',[0 1000 2000],'YTickLabel',{'0','1','2'},'XTickLabel',{'-2','-1','0','1'})






        for k_y = all_data_2(CCI).cycles
            line([k_y k_y],y1,'Color','k','LineStyle','-.','LineWidth',0.001)
        end
cell_title = ['RLS = ',num2str(all_data_2(CCI).age)];
%         
title(cell_title,'FontSize',24);
set(gca,'FontSize',24)
end

