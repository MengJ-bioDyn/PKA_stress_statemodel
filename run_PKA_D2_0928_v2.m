% Exercise 4
% ----------

% This script generates a number of (Gillespie) realizations;
% plots the average and spread of the trajectories (characterized
% by the 2.5 and 97.5 percentiles) against time.
% 

% Ass = d0da*totA/(d0da + k0ka); Bss = d0db*totA/(d0db + k0kb); Css =
% d0dc*totA/(d0dc + k0kc); Sss = k0/d0;
% ASss = k0katotA(d0da + k0ka); Bss = k0kb*totA/(d0db + k0kb); Css =
% k0kc*totA/(d0dc + k0kc);

% Definition of the compete_silencing A
%///////////////////////////////
% Volume
clear all 
% close all
Omega =20;

% 
parain= [1, 1, 5, 6, 5,...
       1,  1, 2.5,...
       .1, 20, 3, 1,...
       .1, 5, 0.2, 1,...
       0.01, 10,  0.1,  10,  0.002, 0.5,  0.5, 8];
   
 parain= [1, 1, 5, 6, 5,...
       1,  1, 2.5,...
       0.1, 20000, 10, 1,...
       .1, .5, 1, 1,...
       0.03, 10,  0.1,  10,  0.002, 2,  0.5, 0.05];
   
% k1 = para(1)*Omega;     
% TotM = para(2)*Omega;   
% km1 = para(3)*Omega;    
% k1r = para(4);   
% km2 = para(5)*Omega; 
% 
% k2 = para(6); 
% TotP = para(7)*Omega; 
% k2r = para(8)/Omega;  
% 
% k3 = para(9)*Omega; 
% a3 = para(10)*Omega;
% b3 = para(11)*Omega;
% d3 = para(12);
% 
% k5 = para(13); 
% a5 = para(14);  
% b5 = para(15); 
% d5 = para(16); 
% 
% a4 = para(17); 
% b4 = para(18)*Omega;  
% k4 = para(19);  
% km4 = para(20);
% 
% d4 = para(21);
% kpd = para(22)*Omega;
% P0 = para(23)*Omega;
% S = para(24)*Omega;


% Get reactions

[ S, h, endSim ] = PKA_D2_0928(parain, Omega);

% Initial condition
% [sR, R, mG, GFP, A, X]


for ini_i =1:1%length(sir2list)
% x0 = round([0 totA 0 totB 0 totC  sir2list(ini_i)]');

% x0 = round(Omega*[ 19 1 1 1 1 0 0]');

x0 = round(Omega*[ 0.9 0.1 0.5 0.1 0.01]'); % for v4b3_0914
% Simulation time
Tmax = 120;
%///////////////////////////////


%% Set up simulation (multiple trajectories)

reset(RandStream.getGlobalStream);

% Number of trajectories to simulate
numReals = 200;

% Initial state and time
t = repmat(0,[1 numReals]);
x = repmat(x0,[1 numReals]);
% load('rdna_p1_ini.mat','P','M','On','Off');
% x = [Off(end,:); On(end,:); M(end,:); P(end,:)];
% clear P M On Off
% Number of steps done
idx = 0;

% Set up arrays to record trajectories
tsampleIdx = ones(1,numReals);
deltaSample = 0.05;
X = nan(size(x,1),ceil(Tmax/deltaSample+1),size(x,2));
Treg = (0:size(X,2)-1)*deltaSample;
T = nan(1,size(X,2),size(x,2));
X(:,1,:) = x;
T(1,1,:) = t;

%% Actual simulation

while any(t<Tmax)
    updateThese = t<Tmax;
    
    if mod(idx,10000)==0
        disp([num2str(sum(updateThese)) ' to update, median(t) = ' num2str(median(t(updateThese)))]);
    end
    idx = idx + 1;
    
    [deltaT, deltaX] = stepGillespie(x(:,updateThese),h,S);
    
    lastx = x;
    lastt = t;

    % Update state and time
    x(:,updateThese) = x(:,updateThese) + deltaX;
    t(updateThese) = t(updateThese) + deltaT;
    
    
    while(true)
        saveThese = t>(1+eps)*tsampleIdx*deltaSample & tsampleIdx+1<=size(X,2);
        if(~any(saveThese))
            break;
        end
        tsampleIdx(saveThese) = tsampleIdx(saveThese)+1;
        for saveThis = find(saveThese)
            X(:,tsampleIdx(saveThis),saveThis) = lastx(:,saveThis);
            T(1,tsampleIdx(saveThis),saveThis) = lastt(saveThis);
        end
    end
    
    % End some simulations
    endThese = endSim(x);
    if any(endThese)==1
        x
    end
    x(:,endThese) = nan;
    t(endThese)=nan;
    
end

% figure;
% Molecule Numbers
%plot(Treg,squeeze(X(1,:,:)));
% Concentrations
%plot(Treg,squeeze(X(1,:,:))/Omega);


%% Plot statistics
colors = {'k','g','m','r','c','b'};
% Comment out for molecule numbers
scaleFactor = 1/Omega;

% X = nan(size(x,1),ceil(Tmax/deltaSample+1),size(x,2));

mean_trace = [];
std_trace =[];
m_lower=[];
m_upper=[];


for idx=1:size(x,1)
%     s_id = soi(idx);
    trace = nanmean(X(idx,:,:),3);
    std_trace = nanstd(X(idx,:,:),1,3);
    lower = prctile(X(idx,:,:),5,3);
    upper = prctile(X(idx,:,:),95,3);
    
    mean_trace = cat(2,mean_trace,trace');
    m_lower = cat(2,m_lower,(trace-std_trace./numReals)');
    m_upper = cat(2,m_upper,(trace+std_trace./numReals)');
    
end

ss = mean_trace(end,:);


% [sR, R, GFP, A, X, D]
for idx =[2 5]
    figure;
    hold on;
    plot(Treg',scaleFactor* mean_trace(:,idx),'linewidth',2); 
    plot(Treg,  scaleFactor*m_lower(:,idx),'r--');
    plot(Treg, scaleFactor* m_upper(:,idx),'r--');
end
    

N= size(X,2);
% si = 70;
% figure; subplot(3,1,1); plot(X(2,:,si)); ylabel('Msn2')
% subplot(3,1,2); plot(X(3,:,si)); ylabel('PKA')
% subplot(3,1,3); plot( X(5,:,si)); ylabel('Damage2'); 


species = cell(1, size(x,1));
for i=1:size(x,1)
    species{i} = squeeze(X(i,:,:));
end



end


%% take threshold
thrshlist=[15 17 20 22 25];
% thrshlist=[ 5 6 7 8 9 10];
thrshN = length(thrshlist);
hist_cell=cell(1, thrshN);

cycT = 2;

D1 = species{5};
GFP= species{5};
GFP_thrsh = zeros(size(GFP));
GFP_thrsh_fold = GFP_thrsh;


for i=1:thrshN
    thrsh_mtx = ones(size(D1))*thrshlist(i);
    delta_D1_thrsh = thrsh_mtx - D1;
    delta_D1_thrsh2 = sign(delta_D1_thrsh);
    
    delta_D1_thrsh2(delta_D1_thrsh2<0)=0;
    
    for tji = 1 : numReals
        t1= find(delta_D1_thrsh2(:,tji)==0,1,'first');
        delta_D1_thrsh2(t1:end,tji)=0;
    end
    GFP_thrsh = GFP.* delta_D1_thrsh2;
    
%     GFP_initial = mean(GFP_thrsh(1:200, :),1);
%     GFP_ini_mtx = repmat(GFP_initial, size(GFP,1) ,1);
%     GFP_fold = GFP_thrsh./GFP_ini_mtx./Omega;
    GFP_fold = GFP_thrsh./Omega;

    
    cut_t = sum(delta_D1_thrsh2,1);
    temp  = Treg(cut_t);
    temp2 = temp.* (temp<Tmax-1);
    hist_cell{i} = temp2(temp2>0);
    figure; subplot(1,2,1); hist(hist_cell{i},[0:cycT:Tmax]);
            subplot(1,2,2); im=imagesc(GFP_fold,[0 2.5]);
%            im.CDataMapping='scaled';colormap('jet')
    [mu sigma] = normfit(hist_cell{i});
    disp('mu sigma');
    [mu sigma]./cycT
end

% good thrsh=22, where[mu sigma]./[24.2 7.5] = [15 15]

%%
% 
 i=4;
 thrshlist(i)
 thrsh_mtx = ones(size(D1))*thrshlist(i);
    delta_D1_thrsh = thrsh_mtx - D1;
    delta_D1_thrsh2 = sign(delta_D1_thrsh);
    
    delta_D1_thrsh2(delta_D1_thrsh2<0)=0;
    
    for tji = 1 : numReals
        t1= find(delta_D1_thrsh2(:,tji)==0,1,'first');
        delta_D1_thrsh2(t1:end,tji)=0;
    end
    GFP_thrsh = GFP.* delta_D1_thrsh2;
    
%     GFP_initial = mean(GFP_thrsh(1:200, :),1);
%     GFP_ini_mtx = repmat(GFP_initial, size(GFP,1) ,1);
%     GFP_fold = GFP_thrsh./GFP_ini_mtx./Omega;
    GFP_fold = GFP_thrsh./Omega;

    
    cut_t = sum(delta_D1_thrsh2,1);
    [cut_t_sort sort_id] = sort(cut_t', 'ascend');
    temp  = Treg(cut_t);
    temp2 = temp.* (temp<Tmax-1);
    hist_cell{i} = temp2(temp2>0);
    GFP_fold2 = GFP_fold(:,sort_id)';
    figure; subplot(1,2,1); hist(hist_cell{i}./cycT,[0:2:50]); xlim([0 50])
            set(gca,'fontsize',14)
            title(['Stress = ',num2str(parain(end))])
%             subplot(1,2,2); im=imagesc(GFP_fold2([11:20,51:60,101:110,151:160,191:200],:),[0 2.5]); xlim([0 12000])
         
          subplot(1,2,2);
%           use_id =[4 6 16:17 20:23 26 27 30:33 35:37 39:41 44:46 48:57 59 61 63:68 70:73 75 77 79 80 84 85 88 91 95 96 99 ]';
          im=imagesc(GFP_fold2,[0 2.5]); xlim([0 Tmax/deltaSample])
%     CC=im.CData;im.CDataMapping='scaled'; colormap('jet')
    xpos =[0 3000 6000]; 
    set(gca, 'XTick', xpos,'XTickLabel',{'0','1000','2000'})
    set(gca,'fontsize',14)
                           
    [mu sigma] = normfit(hist_cell{i});
    title(['mu sigma', num2str([mu sigma]./cycT)])
    
    %%
% %     t_end =3348; trj=68; % S20
% % 
% figure; 
% subplot(2,2,1); plot(1:t_end,GFP_fold(1:t_end, trj),'linewidth',2);
%                 set(gca,'fontsize',14); ylabel('D2 ');  xlabel('Time, min '); ylim([0 3])
%                 title(['Stress = ',num2str(parain(end))])
% subplot(2,2,2); plot(1:t_end,X(2,1:t_end,trj)/Omega,'linewidth',2);
%                 set(gca,'fontsize',14); ylabel('Msn2* ');  xlabel('Time, min '); ylim([0 1.1])
%                 title(['Stress = ',num2str(parain(end))])
% subplot(2,2,3); hist(hist_cell{i}./12,[0:2:50]); xlim([0 50])
%             set(gca,'fontsize',14); ylabel('Counts ');  xlabel('Lifespan ')  
%             title(['mean =  \sigma =', num2str([mu sigma]./12)]);
% % for i=1:2
% %     subplot(2,2,i)
% %      xlabel('Time,min')
% %     xlim([0 12000]); line([0 12000], [2.05 2.05])
% %     xpos =[0 3000 6000 9000 12000 ]; 
% %     set(gca, 'XTick', xpos,'XTickLabel',{'0','1000','2000','3000','4000'})
% %     
% % end
cycT=2;
S=[0.05 0.1 1 2 4 8 12];
for i=1:length(S)
R(i)= load(['run_PKA_D2_v2_S',num2str(S(i)),'_t2.mat'],'hist_cell');
end
X=zeros(199,length(S));
for i=1:length(S)
    X(:,i) = R(i).hist_cell{4}(1:199)';
end
    
figure; boxplot(X/cycT,'notch','on')
set(gca,'fontsize',16)
xlabel('Stress '); ylabel('Lifespan ')

% ta=[R(1).hist_cell{4}]';
% size(ta)
% X(:,1)=[ta;ta(1:97)];
% tb = [R(2).hist_cell{4}]';
% size(tb)
% X(:,2) = [tb; tb(1:46)];
% X(:,2) = [tb; tb(1:36)];
% tb = [R(3).hist_cell{4}]';
% size(tb)
% X(:,3) = [tb; tb(1:76)];
% tb = [R(4).hist_cell{4}]';
% size(tb)
% X(:,4) = [tb; tb(1:21)];
% tb = [R(5).hist_cell{4}]';
% X(:,5) = [tb(1:195)];
% tb = [R(6).hist_cell{4}]';
% X(:,6) = [tb(1:195)];
% tb = [R(7).hist_cell{4}]';
% X(:,7) = [tb; tb(1:95)];
% size(X)
% figure; boxplot(X)
% figure; boxplot(X/12)
% figure; boxplot(X/12,'compact')
% figure; boxplot(X/12,'plotstyle','compact')
% figure; boxplot(X/12,'box')
% figure; boxplot(X/12,'notch','on')
% set(gca,'fontsize',16)
% xlabel('Stress '); ylabel('Lifespan ')