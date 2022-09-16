k=[6 6]; %Points for different Aspect Ratio rectangles
fig=figure;
for q=1:2 %This loop runs through the two points which generate different aspect ratio panels
%Generate a Rectangle with four coordinates
x1=[-1 1];
y1=[-1 1];
[x,y]=meshgrid(x1,y1);
x=reshape(x,[1,4]);
y=reshape(y,[1,4]);
z=zeros(1,4);
% area=area3d(x,y,z);

%This is required to change the limits of the integrations
xlim=sort(unique(x)); 
ylim=sort(unique(y)); 

zp=0:0.1:10;
res_1p=zeros(size(zp,2),1);
res_2p=res_1p;

for i=1:size(zp,2) %This for loop varies the height of the field point
%Field Point
fp=[2 2 zp(i)];

res_1p(i)=gaussq(1,xlim,ylim,fp);

res_2p(i)=gaussq(2,xlim,ylim,fp);
end
subplot(2,1,q)
plot(zp,res_2p-res_1p)
yline(0)
title(sprintf("Aspect Ratio = %d",q))
end
han=axes(fig,'visible','off');
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Difference between 2 point and 1 point Quadrature');
xlabel(han,'Height of the Field Point');
  

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