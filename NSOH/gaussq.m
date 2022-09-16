function res=gaussq(n,xlim,ylim,fp)
syms x

expr=legendreP(n,x); %Generates the Legendre Polynomial
expr_diff=diff(expr,x); %Differentiates it

%Roots will be the same in both X and Y Directions
roots=double(vpasolve(expr==0));

%Weights will be the same in both X and Y direction
weights=2./((1-roots.^2).*double((subs(expr_diff,roots).^2))); 

%Need to do this to perform the double summation
weights_summation=weights*weights';
[rootsx,rootsy]=meshgrid(roots,roots);
rootsx_p=rootsx*diff(xlim)/2+mean(xlim);
rootsy_p=rootsy*diff(ylim)/2+mean(ylim);
%To get the logic of the above equations refer the notes

res=diff(xlim)*0.25*diff(ylim)*sum(sum(weights_summation.*1./sqrt(((rootsx_p-fp(1)).^2)+((rootsy_p-fp(2)).^2)+(fp(3)^2))));
end
% 
% 
% 
% 
% %1 point Gauss Quadrature
% syms x3
% expr1=legendreP(1,x3);
% roots_1=double(vpasolve(expr1==0,x3));
% weights_1=2/((1-roots_1^2)*(subs(diff(expr1,x3),roots_1))^2);
% weights_summation_1=weights_1*weights_1';
% [rootsx_1,rootsy_1]=meshgrid(roots_1,roots_1);
% rootsx_1p=(rootsx_1*diff(xlim)/2)+mean(xlim);
% rootsy_1p=(rootsy_1*diff(ylim)/2)+mean(ylim);
% result_1p=diff(xlim)*diff(ylim)*(1/sqrt(((rootsx_1p)-fp(1))^2+((rootsy_1p)-fp(2))^2+(fp(3))^2));
% 
% 
% %2 Point Gauss Quadrature
% %Legendre Polynomial of 2nd order
% syms x2
% expr=legendreP(2,x2);
% expr_diff=diff(expr,x2);
% roots_2=double(vpasolve(expr==0));%Roots will be the same in both X and Y Directions
% weights_2=2./((1-roots_2.^2).*double((subs(expr_diff,roots_2).^2))); %Weights will be the same in both X and Y direction
% 
% 
% weights_summation=weights_2*weights_2';
% [rootsx_2,rootsy_2]=meshgrid(roots_2,roots_2);
% % diff(xlim)/2,mean(xlim)
% % diff(ylim)/2,mean(ylim)
% rootsx_2p=rootsx_2*diff(xlim)/2+mean(xlim);
% rootsy_2p=rootsy_2*diff(ylim)/2+mean(ylim);
% result_2p=diff(xlim)*0.25*diff(ylim)*sum(sum(weights_summation.*1./sqrt(((rootsx_2p-fp(1)).^2)+((rootsy_2p-fp(2)).^2)+(fp(3)^2))));
% 
