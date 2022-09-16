%Load in N, Length and Area in mm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Nodes of the Truss
nnodes=4; %Number of Nodes
%1st col-Node,2nd Col- X coord, 3rd col - Y coord
nodes=zeros(nnodes,3);
nodes(:,1)=1:nnodes;
nodes(:,2:3)=[750 500;0 500;0 0;0 1000]; %3rd Node is chosen as the origin

%Elements of the truss
nelems=3;
%1st col-Element, 2nd Col-Node1, 3rd Col- Node2, 4th Col- Area, 5th Col - E
%6th col-Length 7th-Orientation of the element
elems=zeros(nelems,5);
elems(:,1)=1:nelems;
elems(:,2:end-1)=[1 3 1250;1 2 1000;1 4 1250];
elems(:,end)=200000*ones(nelems,1); %E of the element
%Computing the length of the elements
elems=[elems zeros(nelems,1) zeros(nelems,1)];
for i=1:nelems
    start=elems(i,2);
    endd=elems(i,3);
    elems(i,end-1)=((nodes(start,2)-nodes(endd,2))^2+(nodes(start,3)-nodes(endd,3))^2)^0.5;
    elems(i,end)=atan((nodes(endd,3)-nodes(start,3))/(nodes(endd,2)-nodes(start,2)))*180/pi;
end


%Loading
%1st col-Node, 2nd col-Fx, 3rd Col- Fy
load=zeros(nnodes,2);
%Load at Node 1
load(:,1)=1:nnodes;
load(1,2)=1000;
load(1,3)=-500;

%Loading Vector
f=zeros(2*nnodes,1);
for i=1:nnodes
    f(2*i-1)=load(i,2);
    f(2*i)=load(i,3);
end

%Forming the global stiffness matrix.
kg=zeros(2*nnodes); %Global Stiffness Matrix
for i=1:nelems
    ael=elems(i,4)*elems(i,5)/elems(i,6); %AE/L
    alpha=elems(i,end); %Angle
    
    
    kl=zeros(4);
    %Local Stiffness Matrix
    kl(1:2:end,1:2:end)=ael*[1 -1;-1 1];
    
    %Transformation Matrix
    t=zeros(4);
    l=cosd(alpha);
    m=sind(alpha);
    t(1:2,1:2)=[l -m;m l];
    t(3:4,3:4)=[l -m;m l];
    
    %Tansformed Stiffness Matrix
    kl_transform=t*kl*t';
    
    %The elements of the above matrix have to be distributed to the Global
    %Stiffness Matrix
    start=elems(i,2);
    endd=elems(i,3);
    
    i1=2*start-1;
    i2=2*start;
    i3=2*endd-1;
    i4=2*endd;
    qw=[i1;i2;i3;i4];
    
    kg(qw,qw)=kg(qw,qw)+kl_transform;
end


%Displacement Boundary Conditions
%Modify the following lines according to the problem
%1st col-Node, 2nd col-X displacement, 3rd col-Y Displacement
%Nodes-2,3,4 have zero displacements in both X and Y direction
dispBC=zeros(nnodes-1,3);
dispBC(:,1)=2:nnodes;

%Penalty Method
kg1=kg;
f1=f;
max_kg=max(max(kg1))*1e4;
for i=1:size(dispBC,1)
        ind=2*dispBC(i,1);
        kg1(ind,ind)=kg1(ind,ind)+max_kg;
        f1(ind)=f1(ind)+max_kg*dispBC(i,3);
        
        kg1(ind-1,ind-1)=kg1(ind-1,ind-1)+max_kg;
        f1(ind-1)=f1(ind-1)+max_kg*dispBC(i,2);
end
U1=kg1\f1;


%Elimination Method
% Therefore the stiffness matrix elements associated (both X and Y) with the three nodes
% have to be eliminated.
for i=size(dispBC,1):-1:1
    if dispBC(i,3)==0
        kg(2*dispBC(i,1),:)=[];
        kg(:,2*dispBC(i,1))=[];
        f(2*dispBC(i,1))=[];
    end
    if dispBC(i,2)==0
        kg(2*dispBC(i,1)-1,:)=[];
        kg(:,2*dispBC(i,1)-1)=[];
        f(2*dispBC(i,1)-1)=[];
    end
end



%Displacement at the the required node
U=kg\f;
fprintf("The X-displacement of the node 1, through elimination method, is %f mm\n",U(1))
fprintf("The Y-displacement of the node 1, through elimination method, is %f mm\n",U(2))
fprintf("The X-displacement of the node 1, through Penalty method, is %f mm\n",U1(1))
fprintf("The Y-displacement of the node 1, through Penalty method, is %f mm\n",U1(2))


    


