function dy = ode_pka_D2_t2(t,y, para)

k1 = para(1);     
TotM = para(2);   
km1 = para(3);    
k1r = para(4);   
km2 = para(5); 

k2 = para(6); 
TotP = para(7); 
k2r = para(8);  

k3 = para(9); 
a3 = para(10);
b3 = para(11);
d3 = para(12);

k5 = para(13); 
a5 = para(14);  
b5 = para(15); 
d5 = para(16); 

a4 = para(17); 
b4 = para(18);  
k4 = para(19);  
km4 = para(20);
d4 = para(21);
kpd = para(22);
P0 = para(23);
S = para(24);

  
dy = zeros(5,1);
%y1,y2,y3, y4
%[M2a P Oxg DD2];
    
    %M2a % I kept the MM kinetics for Msn2, as its activation/deactivation
    % depends on phosphorylation/dephosphorylation
    dy(1) = k1*(TotM - y(1))/(km1 + TotM - y(1)) - k1r*y(2)*y(1)/(km2 + y(1));
    
    %PKA
    dy(2) = k2*(TotP - y(2)) - k2r*y(2)*(S + y(5)) ;

    % Sip18
    dy(3) = k3 + a3*y(1)^4/(b3^4 + y(1)^4)  - d3*y(3); 
    
    % PNC1
    dy(4) = k5 + a5*y(1)^2/(b5^2 + y(1)^2)  - d5*y(4); 
    
    % Cellular damage
    dy(5) = a4/( 1 +  (y(3)/b4)^2) *(  k4*S + kpd*y(2)/(P0 + y(2))) - d4*y(5) ;
   

% end
end


