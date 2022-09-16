function area=calcarea2d(pts1,pts2)
%pts1 and pts2 can be in the combination of XY, YZ and ZX
%First sort them using the centroid
centroid=[mean(pts1),mean(pts2)];

%Now calculate the angle made by each of these points with the centroid
if ~(sum(pts1)==pts1(1)*size(pts1,2)||sum(pts2)==pts1(2)*size(pts1,2)) %If
    angle=zeros(size(pts1));
    for i=1:size(pts1,2)
        a1=pts1(i)-centroid(1);
        a2=pts2(i)-centroid(2);
        if a1>0 && a2>0
            angle(i)=atand(a2/a1);
        elseif a1>0 && a2<0
            angle(i)=360+atand(a2/a1);
        else
            angle(i)=180+atand(a2/a1);
        end
    end
    pts=[pts1;pts2;angle];
    %Sorts these points based on the angle. We get a anti-clockwise sense of
    %the points
    pts=sortrows(pts',3)';
    
    x=pts(1,:);
    y=pts(2,:);
    x_modi=[x,x(1)];
    y_modi=[y,y(1)];
    y_sum=y_modi(1:end-1)+y_modi(2:end);
    x_diff=diff(x_modi);  %Calculates the differences of consecutive elements
    area=abs(0.5*x_diff*y_sum');
else
    area=0;
end
end