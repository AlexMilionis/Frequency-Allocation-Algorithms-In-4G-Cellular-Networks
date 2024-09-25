function coord=subscribers_coordinates(cx,cy,radius,UEspercell)
    
    rng('shuffle')
    corners=zeros(6,2);
    coord = rand(UEspercell,2);    
    h = (sqrt(3)*radius)/2;

    corners(1,1)=cx-radius;      corners(1,2)=cy;
    corners(2,1)=cx-radius/2;    corners(2,2)=cy+h;
    corners(3,1)=cx+radius/2;    corners(3,2)=cy+h;
    corners(4,1)=cx+radius;      corners(4,2)=cy;
    corners(5,1)=cx+radius/2;    corners(5,2)=cy-h;
    corners(6,1)=cx-radius/2;    corners(6,2)=cy-h;



    m=zeros(1,4);     %y=m*x+b
    b=zeros(1,4);

    m(1)=(corners(2,2)-corners(1,2))/(corners(2,1)-corners(1,1));
    b(1)=corners(1,2)-m(1)*corners(1,1);
    m(2)=(corners(6,2)-corners(1,2))/(corners(6,1)-corners(1,1));
    b(2)=corners(1,2)-m(2)*corners(1,1);
    m(3)=(corners(4,2)-corners(5,2))/(corners(4,1)-corners(5,1));
    b(3)=corners(5,2)-m(3)*corners(5,1);
    m(4)=(corners(4,2)-corners(3,2))/(corners(4,1)-corners(3,1));
    b(4)=corners(3,2)-m(4)*corners(3,1);


    for i=1:UEspercell
        coord(i,1)=(2*radius).*coord(i,1)+(cx-radius);   %(b-a).*x + a
        coord(i,2)=(2*h).*coord(i,2)+(cy-h); 
    end


    for i=1:UEspercell
        x0 = coord(i,1);
        y0 = coord(i,2);

        if x0>=cx-radius && x0<cx-radius/2 && y0>cy  && y0>m(1)*x0+b(1)    %&& y0>((2*h)/radius)*x0+2*h
            d1 = (y0-b(1))/m(1);
            d = abs(d1-x0);
            coord(i,1) = x0 + 2*d;
        end
        if x0>=cx-radius && x0<cx-radius/2 && y0<cy &&  y0<m(2)*x0+b(2)  %y0<((-2*h)/radius)*x0-2*h
            d1 = (y0-b(2))/m(2);
            d = abs(d1-x0);
            coord(i,1) = x0 + 2*d;
        end
        if x0>cx+radius/2 && x0<=cx+radius && y0<cy && y0<m(3)*x0+b(3) %y0<((2*h)/radius)*x0-2*h
            d1 = (y0-b(3))/m(3);
            d = abs(d1-x0);
            coord(i,1) = x0 - 2*d;
         end
        if x0>cx+radius/2 && x0<=cx+radius && y0>cy && y0>m(4)*x0+b(4) %y0>(-(2*h)/radius)*x0+2*h
            d1 = (y0-b(4))/m(4);
            d = abs(d1-x0);
            coord(i,1) = x0 - 2*d;
        end
    
    end
    scatter(coord(:,1),coord(:,2),5)
end