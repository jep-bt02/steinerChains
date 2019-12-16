clc
close all

% Complex Analysis Project
% Steiner chains inversion through Möbius transformation in the complex plane
% ---------------------------------------------------------------------------

% This code generates a Steiner chain of "n" circles and uses a Möbius transformation to generate a mapping in the complex plane through a center of inversion and a set of rescaling parameters. Central circle is the unit circle.  

% The arangement of the circles and their radius in the Steiner chain is:
% O---(r_internal)---(r_external)---(r)---(R)-------------

% The math begind the code can be found in these sources:
%   https://ima.org.uk/8218/exploring-steiner-chains-mobius-transformations/
%   https://www.youtube.com/watch?v=sG_6nlMZ8f4

% PARAMETERS ----------------------------------------------------------------
% To change the number of circles in the chain, change the "n" parameter: 

n = 9; % number of circles in the chain
points = 1000;
t = linspace(0,2*pi,points); % angular space
theta = pi/n; % angle between center and tangent point of contact

r_internal = 1; % radius for central circle
r_external = (r_internal*(sin(theta)))/(1-sin(theta)); % radius for circles along the chain 
r = r_internal + r_external; % radius for circle through tangential points
R = r_internal + 2*r_external; % radius for outermost circle

centerOfInversion = -r_external; 



% CIRCLES: Steiner ----------------------------------------------------------
d_centers = zeros(1,n); % initialization of centers array
z_circles = zeros([points n]); % initialization of circles matrix
z0 = exp(i*t); % Unit circle
zR = R*exp(i*t); % External circle
zC = r*exp(i*t); % circle through tangential points

for j = 1:n
    d_centers(j) = centers(theta,j,r); % calcuting 
end

for j = 1:n
    for k = 1:points
        z_circles(k,j) = z0(k)*r_external + d_centers(j);
    end
end

% plots: 
figure; hold on; axis equal; grid on; plot(d_centers,'.'); 
cir1 = plot(z0,'r','LineWidth',2.5); L1 = 'Unit';
cir2 = plot(zR,'b','LineWidth',2.5); L2 = 'Outermost';
cir3 = plot(zC,':','LineWidth',1.5); L3 = 'Chain';
for j = 1:n
    plot(z_circles(1:end,j),'k','LineWidth',1.5)
end
axis([-2.5 3.5 -2.5 2.5]); legend([cir1,cir2,cir3],L1,L2,L3); title("Steiner Chains");



% MAPPING: Steiner Inverted -------------------------------------------------
c1 = 0; c2 = -1; c3 = 0.15; c4 = centerOfInversion; % rescaling parameters

w_centers = zeros(1,n); % initialization of centers array
w_circles = zeros([points n]); % initialization of mapping matrix
w0 = (c1*z0 + c2)./(c3*z0 + c4); % Unit circle mapping
wR = (c1*zR + c2)./(c3*zR + c4); % External circle mapping
wC = (c1*zC + c2)./(c3*zC + c4); % circle through tangential points

for j = 1:n
    for k = 1:points
        w_circles(k,j) = (c1*z_circles(k,j) + c2)./(c3*z_circles(k,j) + c4);
    end
end

for j = 1:n
    w_centers(j) = (c1*d_centers(j) + c2)./(c3*d_centers(j) + c4) ;
end

% plots: 
figure; hold on; axis equal; grid on; plot(w_centers,'.'); 
cir1 = plot(w0,'r','LineWidth',2.5); L1 = 'Unit';
cir2 = plot(wR,'b','LineWidth',2.5); L2 = 'Outermost';
cir3 = plot(wC,':','LineWidth',1.5); L3 = 'Chain';
for j = 1:n
    plot(w_circles(1:end,j),'k','LineWidth',1.5)
end
axis([1 6 -2 2]); legend([cir1,cir2,cir3],L1,L2,L3); title("Mobius transformation");









