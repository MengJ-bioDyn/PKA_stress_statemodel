clear all
% clc
pset= [1, 1,10, 2,0.5, ...
     0.01,0.1, 1.5, 0.2, ... 
    0.005,0.1, 0.2, 0.075, ...
    0.01, 0.2, 0.01, 0.5, ...
    0.2, 0.05, 0.25, 0.5, 0.2, 0.01, 1];

% Omega =1;
% k1 = para(1);     % basal unsilencing
% TotM = para(2);    % basal silencing
% km1 = para(3);      % positive feedback unsilencing
% k1r = para(4);     % 
% km2 = para(5);      % nfb silencing
% k2 = para(6);      % unsilencing ->A
% a2 = para(7);      % A degrading
% b2 = para(8);      % A -> X
% d2 = para(9);     % X degrading
% k3 = para(10);
% a3 = para(11);
% b3 = para(12);
% d3 = para(13);   % mRNA syn
% 
% k4 = para(14);   % GFP syn
% km4 = para(15);   % GFP deg
% d4a = para(16); % kmg for GFP
% d4b = para(17);  % km4
% a4 = para(18);   % syn x
% b4 = para(19);         % deg x
% 
% P0 = para(20);
% g0 = para(21);
% bg = para(22);
% 
% dd = para(23);
% kd = para(24);
% 
% S = para(25);



%y1,y2,y3, y4,y5,y6, y7 
%[sR R GFP A X D O D2];
yini = 1*[0.1 1 0.01 0.01];


T0 = [0:1:500]; %time
options = odeset('RelTol',1e-8, 'AbsTol', 1e-8);
sol0 = ode45(@ode_pka_D2_t1,T0,yini, options,pset);  % solve the ODEs
solution0=deval(sol0,T0);  % evaluate the solutions at timepoints T0
                           % solve the result into matrix solution0

% sol_cmp = ode45(@syn_deg_mono,T0,yini, options,pset_cmp);  % solve the ODE
% solution1=deval(sol_cmp,T0); 
figure;
subplot(3,1,1)
plot(T0,solution0(1,:))
legend('Msn2');

subplot(3,1,2)
plot(T0,solution0(2,:))
legend('PKA');

subplot(3,1,3)
plot(T0,solution0(4,:));
legend('D2');

disp('ss of  D2: ')
[solution0(4,end)]
