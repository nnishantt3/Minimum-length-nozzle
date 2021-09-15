function [v,theta,KL,KR] = calculation(theta_max,theta_O,n)
         nodes = (n*(n+3))/2; % Total number of nodes for all characteristics including wall points
         dt = (theta_max - theta_O)/(n-1) ; % increment of theta angle
         v = zeros(1,nodes);
         theta = zeros(1,nodes);
         KL = zeros(1,nodes);
         KR = zeros(1,nodes);
         
         % Calculation of above parameters for the first reflected
         % characteristic
         for i=1:n
             theta(i) = theta_O + ((i-1)*dt);
             v(i) = theta(i);
             KL(i) = theta(i) - v(i);
             KR(i) = theta(i) + v(i);
         end
         i = n+1 % wall node for first reflected characteristic
         theta(i) = theta(i-1);
         v(i) = v(i-1);
         KL(i) = KL(i-1);
         KR(i)  = KR(i-1);
         
         % calculation of the parameters in each node considering one
         % reflected characteristic at a time
         q = n+2;
         p = 2;
         for k = 1:n-1 % k=1 is the second reflected characteristic and k = (n-1) is the last reflected characteristic
             h = q;
             j = p;
             theta(h) = 0;
             v(h) = KR(j);
             KR(h) = KR(j);
             KL(h) = - v(h);
             j = j+1;
             for i =  h+1:(n+q-p)
                 KR(i) = KR(j);
                 KL(i) = KL(i-1);
                 theta(i) = 0.5 * (KL(i) + KR(i));
                 v(i) = 0.5 * (KR(i) - KL(i));
                 j = j+1;
             end    
             if i == (n+q-p)
                h = i+1;
             else
                h = h+1;
             end    
             theta(h) = theta(h-1);
             v(h) = v(h-1);
             KL(h) = KL(h-1);
             KR(h) = KR(h-1);
             q = h+1;
             p = p+1;
         end
end