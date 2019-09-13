clear all;
close all;

% clc


% pset= [0.05, 1, 1.5, 0.6, 0.75,...
%        0.1,  1, 0.25,...
%        0.002, 0.2, 5, 0.015,...
%        0.002, 0.45, 0.4, 0.1,...
%        0.001, 2.5,  0.05,  1,  0.002, 0.5,  0.5, 1];
   
   pset= [1, 1, 5, 6, 5,...
       1,  1, 2.5,...
       .1, 20, 3, 1,...
       .1, 5, 0.2, 1,...
       0.01, 10,  0.1,  10,  0.002, 0.5,  0.5, 1];


% k1 = para(1);      % msn2 reactions
% TotM = para(2);   
% km1 = para(3);    
% k1r = para(4);   
% km2 = para(5); 
% 
% k2 = para(6);     % PKA reactions
% TotP = para(7); 
% k2r = para(8);  
% 
% k3 = para(9);     % Sip18 reactions
% a3 = para(10);
% b3 = para(11);
% d3 = para(12);
% 
% k5 = para(13);     % PNC1 reactions
% a5 = para(14);  
% b5 = para(15); 
% d5 = para(16); 
% 
% a4 = para(17);     % Damage reactions
% b4 = para(18);  
% k4 = para(19);  
% km4 = para(20);
% d4 = para(21);
% kpd = para(22);
% P0 = para(23);
% S = para(24);

 % MSN2    dy(1) = k1*(TotM - y(1))/(km1 + TotM - y(1)) - k1r*y(2)*y(1)/(km2 + y(1));
 % PKA     dy(2) = k2*(TotP - y(2)) - k2r*y(2)*(S + y(5)) ;
 % Sip18   dy(3) = k3 + a3*y(1)^2/(b3^2 + y(1)^2)  - d3*y(3); 
 % PNC1    dy(4) = k5 + a5*y(1)^2/(b5^2 + y(1)^2)  - d5*y(4); 
 % Damage  dy(5) = a4/( 1 + (y(3)/b4)^2) *(  k4*S + kpd*y(2)/(P0 + y(2) ) ) - d4*y(5) ;

%y1,   y2, y3,   y4,  y5
%[Msn2 PKA Sip18 PNC1 Damage];
yini = [0.1 0.5 0.2 1 0.01];

S =[0.1 1];
D2ss = zeros(1,length(S));

T0 = [0:1:100]; %time
options = odeset('RelTol',1e-8, 'AbsTol', 1e-8);

sol_cell= cell(1,2);

figure(1);
hold off;

% S =[0.2:0.2:1, 2:1:20];
% D2ss = zeros(1,length(S));

T0 = [0:1:300];

for i=1:length(S)
pset(end) = S(i);
sol0 = ode45(@ode_pka_D2_t2,T0,yini, options,pset);  % S=0.1
solution0=deval(sol0,T0);  
sol_cell{i} = solution0;

% pset(end) = 1;
% sol1 = ode45(@ode_pka_D2_t1,T0,yini, options,pset);  % S=1
% solution1=deval(sol1,T0);
 


figure(1);

subplot(3,2,1); 
plot(T0,solution0(1,:));
hold on;
legend('Msn2');

subplot(3,2,2); 
plot(T0,solution0(2,:));
hold on;
legend('PKA');


subplot(3,2,3); 
plot(T0,solution0(3,:));
hold on;
legend('Sip18');

subplot(3,2,4); 
plot(T0,solution0(4,:));
hold on;
legend('Pnc1');

subplot(3,2,5);
plot(T0,solution0(5,:));
hold on;
legend('Damage');

pause;

% D2ss(i)=[solution0(4,end)];

% [solution0(:,end)]
end

% figure;
% plot(S, D2ss,'linewidth',2);
% set(gca,'fontsize',14);
% xlabel('Stress '); ylabel('D2 steady state ')

% subplot(3,2,1); 
% legend('Msn2, S=0.1', 'Msn2, S=1');
% subplot(3,2,2); 
% legend('PKA, S=0.1', 'PKA, S=1');
% subplot(3,2,3); 
% legend('Spi18, S=0.1', 'Spi18, S=1');
% subplot(3,2,4); 
% legend('Pnc1, S=0.1', 'Pnc1, S=1');
% subplot(3,2,5);
% legend('Damage');



% get scale factors for msn2, Sip18, PNC1
% use 10% of time, T=30
M_scl = sol_cell{1}(1,31); 
Sip18_scl = sol_cell{1}(3,31); 
PNC1_scl = sol_cell{1}(4,31); 
% 
figure(3);
hold off;
subplot(3,1,1); plot(T0, sol_cell{1}(1,:)./M_scl,'r--', T0, sol_cell{2}(1,:)./M_scl,'r-');
set(gca,'fontsize',14);
xlabel('Stress '); ylabel('Msn2 ')
subplot(3,1,2); plot(T0, sol_cell{1}(3,:)./Sip18_scl,'c--', T0, sol_cell{2}(3,:)./Sip18_scl,'c-');
set(gca,'fontsize',14);
xlabel('Stress '); ylabel('Sip18 ')
subplot(3,1,3); plot(T0, sol_cell{1}(4,:)./PNC1_scl,'g--', T0, sol_cell{2}(4,:)./PNC1_scl,'g-');
set(gca,'fontsize',14);
xlabel('Stress '); ylabel('Pnc1')
