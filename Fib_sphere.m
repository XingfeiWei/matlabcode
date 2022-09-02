num_points=400;
    gr=(sqrt(5.0) + 1.0) / 2.0; % golden ratio = 1.6180339887498948482
    ga=(2.0 - gr) * (2.0*pi);  % golden angle = 2.39996322972865332

for i=1:num_points
        lat = asin(-1.0 + 2.0 * double(i) / (num_points+1));
        lon = ga * i;

        x = cos(lon)*cos(lat);
        y = sin(lon)*cos(lat);
        z = sin(lat);

        result(i,1)=x;result(i,2)=y;result(i,3)=z;
end