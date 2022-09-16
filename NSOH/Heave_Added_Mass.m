l=6;
b=3;
h=2;

%There will be 6 faces: lb, bh, lh;
%No. of faces on each of them depend on the area of each face. we divide
%each of them into unit area panels:
%There will be two sets of each of these panels
n_xy=l*b;
n_yz=b*h;
n_zx=l*h;
n_panels=(n_xy+n_yz+n_zx)*2;

%1st-X,2nd-Y,3rd-Z,4th-nz
cg=zeros((n_xy+n_yz+n_zx)*2,3);
n=zeros((n_xy+n_yz+n_zx)*2,3);

%Lets take the lower-left corner on the face away from us to be the origin.
points_l=0.5:1:l;
points_b=0.5:1:b;
points_h=0.5:1:h;

%The third dimenssion of the 
%1-X Coordinate, 2- Y Coordinate, 3- Z Coordinate
[cg_xy_x,cg_xy_y]=meshgrid(points_l,points_b);
for i=1:2*n_xy
    if i<=n_xy
        cg(i,:)=[cg_xy_x(i) cg_xy_y(i) 0];
        n(i,:)=[0 0 -1];
    else
        cg(i,:)=[cg_xy_x(i-n_xy) cg_xy_y(i-n_xy) h];
        n(i,:)=[0 0 +1];
    end
end

[cg_yz_y,cg_yz_z]=meshgrid(points_b,points_h);
for i=1:2*n_yz
    if i<=n_yz
        cg(i+2*n_xy,:)=[0 cg_yz_y(i) cg_yz_z(i)];
        n(i+2*n_xy,:)=[-1 0 0];
    else
        cg(i+2*n_xy,:)=[l cg_yz_y(i-n_yz) cg_yz_z(i-n_yz)];
        n(i+2*n_xy,:)=[+1 0 0];
    end
end

[cg_xz_x,cg_xz_z]=meshgrid(points_l,points_h);
for i=1:2*n_zx
    if i<=n_zx
        cg(i+2*n_xy+2*n_yz,:)=[cg_xz_x(i) 0 cg_xz_z(i)];
        n(i+2*n_xy+2*n_yz,:)=[0 -1 0];
    else
        cg(i+2*n_xy+2*n_yz,:)=[cg_xz_x(i-n_zx) b cg_xz_z(i-n_zx)];
        n(i+2*n_xy+2*n_yz,:)=[0 +1 0];
    end
end

A=2*pi*eye(size(cg,1));
C=zeros(size(cg,1),1);
   
for i=1:n_panels
    p=[cg(i,1) cg(i,2) cg(i,3)];
    C_sum=0;
    for j=1:n_panels
        if i~=j
            q=[cg(j,1) cg(j,2) cg(j,3)];
            r=p-q;
            A(i,j)=sum(r.*n(j,:))/(norm(r))^3;
            C_sum=C_sum+n(j,3)/norm(r);
        end
    end
    C(i)=C_sum;
end
phi=A\C;
addedmass=sum(phi(1:2*n_xy).*n(1:2*n_xy,3))*1.025;
fprintf("The dimensions of the cuboid are %d m, %d m, %d m\n",l,b,h)
fprintf("The Heave added mass of the cuboid is %f T\n",addedmass)






