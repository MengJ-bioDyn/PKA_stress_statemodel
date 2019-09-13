function [ S, h, endSim ] = PKA_D2_0928(para, Omega)
% model PKA regulated D2


k1 = para(1)*Omega;     
TotM = para(2)*Omega;   
km1 = para(3)*Omega;    
k1r = para(4);   
km2 = para(5)*Omega; 

k2 = para(6); 
TotP = para(7)*Omega; 
k2r = para(8)/Omega;  

k3 = para(9)*Omega; 
a3 = para(10)*Omega;
b3 = para(11)*Omega;
d3 = para(12);

k5 = para(13); 
a5 = para(14);  
b5 = para(15); 
d5 = para(16); 

a4 = para(17); 
b4 = para(18)*Omega;  
k4 = para(19);  
km4 = para(20);

d4 = para(21);
kpd = para(22)*Omega;
P0 = para(23)*Omega;
Sig = para(24)*Omega;



% Rate constants [m: total reactions]
c = @(x)([k1*ones(1,size(x,2))
          k1r*ones(1,size(x,2))
          k2*ones(1,size(x,2))
          k2r*ones(1,size(x,2))
          k3*ones(1,size(x,2))
          a3*ones(1,size(x,2))
          d3*ones(1,size(x,2))
          k4*ones(1,size(x,2))
          kpd*ones(1,size(x,2))
          d4*ones(1,size(x,2)) ]);
         
% Stoichiometry matrix [n by m] row: species, col: total reactions 12
% R1: k1*S*(TotM - y(2))/(km1 + TotM - y(2))
% R2: - k1r*y(3)*y(2)/(km2 + y(2))

% P  k2*(TotP - y(3)) - k2r*y(3)*(S + y(5)) ;
% R3: k2*(TotP - y(3))
% R4: - k2r*y(3)*(S + y(5))
% 
% Ox  % k3 + a3*y(1)^2/(b3^2 + y(1)^2)  - d3*y(3); 

% R5: k3
% R6: a3*y(2)^2/(b3^2 + y(2)^2)
% R7: - d3*y(3)

% D2  % a4/( 1 + (y(3)/b4)^2) *(  k4*S + kpd*y(2)/(P0 + y(2))) - d4*y(5)

% R8: a4/( 1 + (y(4)/b4)^2) * k4*S 
% R9: a4/( 1 + (y(4)/b4)^2) * kpd*y(3)/(P0 + y(3))
% R10: - d4*y(5)
% 

   %  R1 R2  R3  R4 R5  R6 R7  R8  R9 R10 R11 R12 R13 R14 R15 R16 R17 R18 R19 R20 R21 R22  R23 R24
S = [ -1   1  0   0  0    0   0   0  0    0                 % 1 inactive Msn2
       1  -1  0   0  0    0   0   0  0    0                  % 2 active Msn2
       0  0   1   -1 0    0   0   0  0    0                  % 3 P
       0  0   0   0  1    1  -1   0  0    0                  % 4 Ox
       0  0   0   0  0    0   0   1  1   -1  ];            % 5 D2
      
  

% c = @(x)([k1*ones(1,size(x,2))
%           k1r*ones(1,size(x,2))
%           k2*ones(1,size(x,2))
%           k2r*ones(1,size(x,2))
%           k3*ones(1,size(x,2))
%           a3*ones(1,size(x,2))
%           d3*ones(1,size(x,2))
%           k4*ones(1,size(x,2))
%           kpd*ones(1,size(x,2))
%           d4*ones(1,size(x,2)) ]);

      

% Stoichiometry matrix [n by m] row: species, col: total reactions 12
% R1: k1*S*y(1)/(km1 + y(1))
% R2: - k1r*y(3)*y(2)/(km2 + y(2))

% P  
% R3: k2*(TotP - y(3))
% R4: - k2r*y(3)*(S + y(5))
 
% Ox  
% R5: k3
% R6: a3*y(2)^2/(b3^2 + y(2)^2)
% R7: - d3*y(4)

% D2  
% R8: a4/( 1 + (y(4)/b4)^2) * k4*S 
% R9: a4/( 1 + (y(4)/b4)^2) * kpd*y(3)/(P0 + y(3))
% R10: - d4*y(5)


h = @(x)(c(x).*[ x(1,:)./(km1 + x(1,:))        % S-> Msn2*
                 x(3,:).*x(2,:)./(km2 + x(2,:))  % P-> Msn2
                 % PKA
                 (TotP*ones(1,size(x,2)) - x(3,:))
                 x(3,:).*(Sig + x(5,:))
                 % Sip18
                 ones(1,size(x,2))
                 x(2,:).^4./(b3^4 + x(2,:).^4)
                 x(4,:)
                 % D2
                 a4*Sig*ones(1,size(x,2))./( 1 + (x(4,:)/b4).^2)
                 a4*ones(1,size(x,2))./( 1 + (x(4,:)/b4).^2).*(x(3,:)./(x(3,:)+P0) ) 
                 x(5,:) ]);
             


% Never end the simulation
endSim = @(x)(false(1,size(x,2)));


end

