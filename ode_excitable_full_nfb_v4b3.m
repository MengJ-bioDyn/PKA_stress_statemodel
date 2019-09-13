function dy = ode_excitable_full_nfb_v4b3(t,y, para)

Omega =1;
kon = para(1);     % basal unsilencing
koff = para(2);    % basal silencing
kp = para(3);      % positive feedback unsilencing
km1 = para(4)*Omega*Omega;     % 
kn = para(5);      % nfb silencing
k2 = para(6)*Omega;      % unsilencing ->A
d2 = para(7);      % A degrading
k3 = para(8)*Omega;      % A -> X
d3a = para(9);     % X degrading
TotR = para(10)*Omega;
km3 = para(11)*Omega*Omega;
km2 = para(12)*Omega*Omega;
% k0 = para(13);   % mRNA syn
% dm = para(14);   % mRNA deg
kg = para(14)*Omega;   % GFP syn
dg = para(15);   % GFP deg
kmg = para(16)*Omega*Omega; % kmg for GFP
km4 = para(17)*Omega;  % km4
sx = para(18)*Omega;   % syn x
dx = para(19);         % deg x

sd = para(20)*Omega;
dd = para(21);
kmd = para(22)*Omega*Omega;

ko = para(23)*Omega;
ko2 = para(24)*Omega*Omega*Omega;
kor = para(25);
kmo = para(26)*Omega*Omega;

sd = para(27)*Omega;
dd = para(28);
kmd = para(29)*Omega*Omega;

  
dy = zeros(6,1);
%y1,y2,y3, y4,y5,y6, y7, y8
%[sR R GFP A X D O, D2];

% if mod(y(4),cycT) > dilut_t
%     dy(1) =  kon*(Tot - y(1)) + kp*y(1)*y(1)*(Tot-y(1))/(km1 + Tot- y(1)) - koff*y(1) - kn*y(3)^2*y(1) ; % rDNA
%     dy(2) = k2*y(1) - d2a*y(2) - d2*y(2);    % A
%     dy(3) = k3*y(2) - d3*y(3);               % X
%     dy(4) =  alpha;                          % Time 
%     dy(5) = k4*y(1) - d4a*y(5) - d4b*y(5);
% else
    dy(1) = - kon*y(1) - kp*y(2)^2*y(1)/(km1 + y(2)^2) + koff*y(2) + kn*y(5)*y(2)/(y(2)+y(5));
    dy(2) =  kon*y(1) + kp*y(2)^2*y(1)/(km1 + y(2)^2) - koff*y(2) - kn*y(5)*y(2)/(y(2)+y(5));

    % GFP reporter
    dy(3) = kg*y(2)^2/(kmg+y(2)^2) - dg*y(3); 
    
    % A
    dy(4) = k2*y(2)^2/(km2+ y(2)^2) - d2*y(4);
    % X
    dy(5) = sx - dx*y(5) + k3*y(4)^2/(km3+y(4))  - d3a*y(5);
    % D1
    dy(6) = sd*y(2)^2/(kmd+y(2)^2) - dd*y(6); 
    
    % O
    dy(7) = ko + ko2/(kmo + y(5)^2) - kor*y(7);
    % D2
    dy(8) = sd2*y(7)^2/(kmd2+y(7)^2) - dd2*y(8); 

% end
end


