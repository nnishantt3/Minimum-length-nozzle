% Input Parameters
Me = input('Enter exit Mach number');
TR = input('Enter Nozzle Throat radius ')
g = input('Enter value for ratio of specific heats')
n = input('Enter number of characteristic lines')

% Calculation of maximum wall angle and dividing this angle into number of
% elementary angles
theta_max = (0.5 * Prandtl_Meyer(Me,g)) * (180/pi); % maximum wall angle
disp(theta_max)
theta_O = (theta_max)/n; %first flow deflection angle

% Importing values of v(Prandtl Meyer functional value),theta(flow deflection angle for each characteristic),KL(equation
% for left running characteristic),KR(equation for right running
% characteristic) at all grid points from the function 'calculation'
[v,theta,KL,KR] = calculation(theta_max,theta_O,n);
nodes = (n*(n+3))/2;
x = zeros(1,nodes);
y = zeros(1,nodes);
M = zeros(1,nodes);
mu = zeros(1,nodes);
for i = 1:nodes
    M(i) = Inv_Prandtl_Meyer(v(i)); % calculation of Mach number at each node using inverse Prandtl Meyer function
end
for i = 1:nodes
    mu(i) = meu(M(i)); % calculation of mu value at each node
end

% Calculation and plotting of the grid points for the first reflected characteristic
for i = 1:(n+1)
    if i == 1 % Point on centreline
       x(i) = -(TR)/(tand(theta(i)-mu(i)));
       y(i) = 0;
       plot([0 x(i)],[TR y(i)],'b')
       hold on
    else if i == n+1 % point on wall
            m1 = tand(0.5 * (theta(i) + theta(i-1) + mu(i) + mu(i-1)));
            m2 = tand(0.5 * (theta_max+theta(i)));
            x(i) = ((-TR + y(i-1) - (m1*x(i-1)))/(m2-m1));
            y(i) = (m2*x(i)) + TR;
            plot([0 x(i)],[TR y(i)],'k')
            hold on
            plot([x(i-1) x(i)],[y(i-1) y(i)],'b')
            hold on
    else 
        a1 = tand(0.5 * (theta(i) + theta(i-1) + mu(i) + mu(i-1)));
        a2 = tand(theta(i) - mu(i));
        x(i) = (y(i-1) - (a1 * x(i-1))- TR)/(a2 - a1);
        y(i) = TR + (a2*x(i));
        plot([0 x(i)],[TR y(i)],'b')
        hold on
        plot([x(i-1) x(i)],[y(i-1) y(i)],'b')
        hold on 
    end
    end
    i = i+1; % after last iteration i takes value of (n+2) grid point
end

%% Calculation of the other grid points including wall for each characteristic one by one
h = i;% h = n+2
k = 0
for j=1:(n-1) % j=1 is the second reflected characteristic and j=n-1 is the last reflected characteristic
    for i = h:(h+n-1-k)
        if i == h  % centreline
            m2 = tand(0.5 * (theta(i) + theta(i-n+k) - mu(i) - mu(i-n+k)));
            x(i) = ((m2 * x(i-n+k)) - y(i-n+k))/m2 ;
            y(i) = 0;
            plot([x(i-n+k) x(i)],[y(i-n+k) y(i)],'b')
            hold on
        else if i == h+n-1-k % wall
                m1 = tand(0.5 * (theta(i) + theta(i-1) + mu(i) + mu(i-1)));
                m2 = tand(0.5 * (theta(i) + theta(i-n+k)));
                x(i) = ((m2*x(i-n+k)) - (m1*x(i-1)) + y(i-1) - y(i-n+k))/(m2-m1);
                y(i) = (m2*(x(i) - x(i-n+k))) + y(i-n+k);
                plot([x(i-n+k) x(i)],[y(i-n+k) y(i)],'k')
                hold on
                plot([x(i-1) x(i)],[y(i-1) y(i)],'b')
                hold on
        else % points between centreline and wall
            m1 = tand(0.5 * (theta(i) + theta(i-1) + mu(i) + mu(i-1)));
            m2 = tand(0.5 * (theta(i) + theta(i-n+k) - mu(i) - mu(i-n+k)));
            x(i) = ((m2 * x(i-n+k)) - (m1 * x(i-1)) + y(i-1) - y(i-n+k))/(m2-m1);
            y(i) = y(i-1) + m1*(x(i) - x(i-1));
            plot([x(i-n+k) x(i)],[y(i-n+k) y(i)],'b')
            hold on
            plot([x(i-1) x(i)],[y(i-1) y(i)],'b')
            hold on
        end
        end
        i = i+1;
    end
    h = i;
    k=k+1;
end

xlabel('Nozzle axis (mm)')
ylabel('Radius (mm)')
title('Contour for minimum length nozzle')
axis equal

% Writing coordinates into a text file
%fid = fopen('MLN_Points_x.txt','w');
%fprintf(fid,'%6.2f \n',x);
%fclose(fid);

%fid = fopen('MLN_Points_y.txt','w');
%fprintf(fid,'%6.2f \n',y);
%fclose(fid);

%fid = fopen('MLN_Points_-y.txt','w');
%fprintf(fid,'%6.2f \n',-y);
%fclose(fid);


