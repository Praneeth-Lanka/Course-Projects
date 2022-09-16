%Load in N, Length and Area in mm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a_load=1; %0 for no axial stiffness, 1 for considering axial stiffness


%Nodes of the Beams
nnodes=3; %Number of Nodes
%1st col-Node,2nd Col- X coord, 3rd col - Y coord
nodes=zeros(nnodes,3);
nodes(:,1)=1:nnodes;
nodes(:,2:3)=[0 0;0 3000;7000 4000];%1st Node is chosen as the origin

%Elements of the truss
nelems=2;
%1st col-Element, 2nd Col-Node1, 3rd Col- Node2, 4th Col- Depth 5th-Thickness, 6th Col - E
%7th col-Length 8th-Orientation of the element
elems=zeros(nelems,6);
elems(:,1)=1:nelems;
elems(:,2:3)=[1 2;2 3];
elems(:,4:5)=[45 200;45 200];  
elems(:,6)=200000*ones(nelems,1); %E of the element
%Computing the length of the elements
elems=[elems zeros(nelems,1) zeros(nelems,1)];
for i=1:nelems
    start=elems(i,2);
    endd=elems(i,3);
    elems(i,end-1)=((nodes(start,2)-nodes(endd,2))^2+(nodes(start,3)-nodes(endd,3))^2)^0.5;
    elems(i,end)=atan((nodes(endd,3)-nodes(start,3))/(nodes(endd,2)-nodes(start,2)))*180/pi;
end

%Loading
%1st col-Node, 2nd col-Fx, 3rd Col- Fy 4th Col - Moment
load=zeros(nnodes,4);
%Load at Node 1
load(:,1)=1:nnodes;
load(1,4)=-100e6;
load(2,2)=150e3;
load(2,3)=-150e3;


%Loading Vector
f=zeros(3*nnodes,1);
for i=1:nnodes
    f(3*i-2)=load(i,2);
    f(3*i-1)=load(i,3);
    f(3*i)=load(i,4);
end

%Forming the global stiffness matrix.
kg=zeros(3*nnodes); %Global Stiffness Matrix
for i=1:nelems
    area=elems(i,4)*elems(i,5); %Cross-Sectional Area of the element
    mi=0.25e9; %Moment of Inertia of element
    L=elems(i,7); %Length of the element
    
    ax_sfness=area*elems(i,6)/L; %AE/L
    flex_sfness=elems(i,6)*mi/L^3; %EI/L^3
    
    alpha=elems(i,end); %Orientation of the beam element
    
    kl=zeros(6);
    %Local Stiffness Matrix
    if a_load==1
        kl(1:3:end,1:3:end)=ax_sfness*[1 -1;-1 1];
    end
    kl(2:3,2:3)=flex_sfness*[12 6*L;6*L 4*L^2];
    kl(5:6,5:6)=flex_sfness*[12 -6*L;-6*L 4*L^2];
    kl(5:6,2:3)=flex_sfness*[-12 -6*L;6*L 2*L^2];
    kl(2:3,5:6)=flex_sfness*[-12 6*L;-6*L 2*L^2];
    
    %Transformation Matrix
    t=zeros(6);
    l=cosd(alpha);
    m=sind(alpha);
    t(1:2,1:2)=[l -m;m l];
    t(3,3)=1;
    t(4:5,4:5)=[l -m;m l];
    t(end,end)=1;
    
    %Tansformed Stiffness Matrix
    kl_transform=t*kl*t';
    
    %The elements of the above matrix have to be distributed to the Global
    %Stiffness Matrix
    start=elems(i,2);
    endd=elems(i,3);
    
    i1=3*start-2;
    i2=3*start-1;
    i3=3*start;
    i4=3*endd-2;
    i5=3*endd-1;
    i6=3*endd;
    qw=[i1;i2;i3;i4;i5;i6];
    
    kg(qw,qw)=kg(qw,qw)+kl_transform;
    
end

%Displacement Boundary Conditions
%Modify the following lines according to the problem
%1st col-Node, 2nd col-X displacement, 3rd col-Y Displacement,
% 4th col-Slope
%Node-1 has zero displacement and slope
%Node-3 has zero slope
dispBC=zeros(2,4);
dispBC(:,1)=[1,3];
%NaN is just a way to distinguish between the zero displacement boundary
%condition and other zeros of the dispBC matrix
dispBC(1,2:3)=NaN;
dispBC(2,2:3)=NaN;


%Non-zero Displacement indices
nzero_indices=1:nnodes*3;
for i=size(dispBC,1):-1:1
    if isnan(dispBC(i,4))
        nzero_indices(3*dispBC(i,1))=[];
    end
    if isnan(dispBC(i,3))
        nzero_indices(3*dispBC(i,1)-1)=[];
    end
    if isnan(dispBC(i,2))
        nzero_indices(3*dispBC(i,1)-2)=[];
    end
end
        
flag=0;
%Removes the elements of the stiffness matrix associated with axial
%displacement
if a_load==0
    for i=nnodes:-1:1
        kg(3*i-2,:)=[];
        kg(:,3*i-2)=[];
        f(3*i-2)=[];
    end
    for i=size(dispBC,1):-1:1
        %Slope
        if isnan(dispBC(i,4))
            kg(2*dispBC(i,1),:)=[];
            kg(:,2*dispBC(i,1))=[];
            f(2*dispBC(i,1))=[];
        end
        %Vertical Displacement
        if isnan(dispBC(i,3))
            kg(2*dispBC(i,1)-1,:)=[];
            kg(:,2*dispBC(i,1)-1)=[];
            f(2*dispBC(i,1)-1)=[];
        end
    end
end

if a_load~=0
    for i=size(dispBC,1):-1:1
        %Slope
        if isnan(dispBC(i,4))
            kg(3*dispBC(i,1),:)=[];
            kg(:,3*dispBC(i,1))=[];
            f(3*dispBC(i,1))=[];
        end
        %Vertical Displacement
        if isnan(dispBC(i,3))
            kg(3*dispBC(i,1)-1,:)=[];       
            kg(:,3*dispBC(i,1)-1)=[];
            f(3*dispBC(i,1)-1)=[];
        end
        %Horizontal Displacement
        if isnan(dispBC(i,2))
            kg(3*dispBC(i,1)-2,:)=[];
            kg(:,3*dispBC(i,1)-2)=[];
            f(3*dispBC(i,1)-2)=[];
        end
    end
end
U=kg\f;

for i=1:size(nzero_indices,2)
    a=ceil(nzero_indices(i)/3);
    xdisp=3*a-2;
    ydisp=3*a-1;
    theta=3*a;
    if nzero_indices(i)==xdisp
        fprintf("The Horizontal Displacement at Node-%d is %f mm\n",a,U(i))
    elseif nzero_indices(i)==ydisp
        fprintf("The Vertical Displacement at Node-%d is %f mm\n",a,U(i))
    elseif nzero_indices(i)==theta
        fprintf("The Slope at Node-%d is %f rad\n",a,U(i))
    end
end

