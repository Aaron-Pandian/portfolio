function inCircle = insideCircle(points,circle)
    [rows,~] = size(points);
    inCircle = [];
    i = 1;
    while i <= rows
        x = points(i,1);
        y = points(i,2);
        if ((x - circle(1))*(x-circle(1)) + (y-circle(2))*(y-circle(2)) < circle(3)*circle(3))
            inCircle = [inCircle,1];
        else
            inCircle = [inCircle,0];
        end
        i = i+1;
    end
end

