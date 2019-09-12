%% Plot two colors on the same figure,can also use for single color plot
%% the age sort and death type selection has been commented off
%%
function plot_mc
red = [1,0,0];
green = [0,1,0];

all_data = input('Please type in the data set you want to plot,e.g. ''WT_20151112'': ');%Get the name of data
Flu = input('Please type in the fluorescence channel you want to plot, e.g. [1 2] :');%Get the fluorescence channel number
% Label = input('Please type in the label for the fluorescence channel, e.g. ''NYN'' :');
DT = input('Please type in what death type you want to plot, eg.. 1 : ');
eval(['load ' all_data ';']);%load data
eval(['all_data =' all_data ';']);%Rename the data want to plot to all_data
all_data = nestedSortStruct(all_data,'age');%Sort the data according to their lifespan
all_data = all_data([all_data.Death_type] == DT);%Plot only death 1 and death 2



% Get the max of axis
max_x = zeros(1,length(all_data));%Create a empty array to store max of x axis of each cell
max_y = zeros(2,length(all_data));%Create a empty array to store max of y axis of each cell

for i_xy = 1:length(all_data)
% %         all_data(i_xy).id%for debug
% % 
% %     all_data(i_xy).cycle%for debug
    
    max_x(i_xy) = all_data(i_xy).cycle(end);%Get the last lifespan of the cell
    max_y(1,i_xy) = max(all_data(i_xy).traj_normalized(all_data(i_xy).cycle(5):all_data(i_xy).cycle(end),:,Flu(1)));%Get the max flu during the current lifespan
    max_y(2,i_xy) = max(all_data(i_xy).traj_normalized(all_data(i_xy).cycle(5):all_data(i_xy).cycle(end),:,Flu(2)));
end
max_y1 = max(max_y(1,:));
max_y2 = max(max_y(2,:));
max_x = max(max_x);


% Specify how many each want for each subplot: Image Per Subplot
IPS = input('please specifiy how many subplot per figure: pick 4,9,16,25 or 36 ');
% Get the column and/or row number of each subplot: Column And Row
CNR = sqrt(IPS);
% Get the number of subplot is need: Number of Subplot
NOS = ceil(length(all_data)/IPS);


% Initiating subplot
% For each figure
for FN = 1:NOS
    figure
    %Determin the cell index on the subplot
    if length(all_data) >= IPS*FN
        CI = 1:IPS;
    elseif length(all_data) < IPS*FN
        CI = 1:(length(all_data)-IPS*(FN-1));
    end
    
    for i_CI = CI
        subplot(CNR,CNR,i_CI);
        %Get the Current Cell Index
        CCI = (FN-1)*IPS + i_CI;
        %Get the end frame of cell lifespan
        life_end = all_data(CCI).cycle(end);
        %Get the start frame of cell lifespan
        life_start = all_data(CCI).cycle(3);
        %Get the frames of lifespan for index fluorescence
        FLS = life_start:life_end;
        %Get the cell cycle info
        cycles = all_data(CCI).cycle(3:end);
        %Remove the zeros
        cycles = cycles(cycles>0);
  
        %Extract the curr_traj info and smooth it within 6 window
        curr_trace_Flu = all_data(CCI).traj_normalized(FLS,:,Flu(1));
        curr_trace_End = all_data(CCI).traj_normalized(FLS,:,Flu(2));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% To convert the unit of x axis from frames to minius %%%%%%%
        date = all_data(CCI).Date;
        date = num2str(date);
        if strcmp(date,'20151112')
            FLS = 6*FLS;
            curr_trace_Flu = smooth(curr_trace_Flu,30);
            curr_trace_End = smooth(curr_trace_End,30);
            cycles2 = 6*cycles;
            
        elseif strcmp(date,'20160309')
            FLS = 6*FLS;
            curr_trace_Flu = smooth(curr_trace_Flu,30);
            curr_trace_End = smooth(curr_trace_End,30);
            cycles2 = 6*cycles;
        else
            FLS = 15*FLS;
            curr_trace_Flu = smooth(curr_trace_Flu,12);
            curr_trace_End = smooth(curr_trace_End,12);
            cycles2 = 15*cycles;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
        %Plot Fluorescence channel and the Nuclear marker Channel
        [hAx,hl1,hl2] = plotyy(FLS,curr_trace_Flu,FLS,curr_trace_End);
        
        box off
        
        %Set up the properties of 1st axis
%         hAx(1).YLim = [0 2];        
        hAx(1).YTickMode = 'auto';
%         hAx(1).XLim = [0 15*max_x];%% use 6 min interval or 15 min interal depend on the data
        hl1.LineWidth = 1.5;
        hAx(1).YColor = 'blue';
        hl1.Color = 'blue';
        C1 = num2str(Flu(1));
        C1 = ['Channel ',C1];
        ylabel(hAx(1),C1);
        
        
        %Set up the properties of 2nd axis
%         hAx(2).YLim = [0 5];
        hAx(2).YTickMode = 'auto';
%         hAx(2).XLim = [0 15*max_x];
        hl2.LineWidth = 1.5;
        hAx(2).YColor = 'red';
        hl2.Color = 'red';
        C2 = num2str(Flu(2));
        C2 = ['Channel ',C2];
        ylabel(hAx(2),C2);

        xlabel(hAx(1),'Time (min)');
       
        cell_title = [num2str(all_data(CCI).age) ',',all_data(CCI).id, ',',num2str(all_data(CCI).Date),',','DT', num2str(all_data(CCI).Death_type)];
        
        title(cell_title);
        
        %Add the cell cycle line
        y1 = get(gca,'ylim');
        for k_y = cycles2
            line([k_y k_y],y1,'Color','k','LineStyle','--','LineWidth',1.5)
        end
        set(gca,'FontSize',12)
        
    end
    
% % % % % % %     [ax,h3]=suplabel('NTS1-pGPD-GFP' ,'t'); 
% % % % % % %     set(h3,'FontSize',20)
   
end

