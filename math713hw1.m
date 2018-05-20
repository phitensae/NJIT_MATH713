%% homework for math 713 - gibbs phenomenon part 1; - demonstrating that the order of convergence is 1/2 in L^2
clear all
close all
clc
clf

%h = 1;
dx = [.5 .25 .125 .1 .01 .005];
N = length(dx); 
xmax = 10; 
clf 
E = zeros(size(dx)); %this will hold the error 
square = @(x) (abs(x) <=1); 
locations = zeros(size(dx)); %this will hold the locations of where the maximum overshoot occurs 
overshoots = zeros(size(dx)); %will contain the information on how much overshoot there was 

for kk = 1:N
    h = dx(kk); %set the stepsize 
    x = -xmax:h:xmax; %computational grid 
    xx = -xmax-h/20:h/10:xmax+h/20; %plotting grid
    v = (abs(x) <=1); %v is the values for a square wave from -1 to 1
    p = zeros(size(xx)); 
    for ii = 1:length(x)
        p = p + v(ii) * h/(2*pi) .* cos((xx - x(ii))/2) .*  sin(pi*(xx - x(ii))/h)./sin((xx-x(ii))/2) ;
    end
    subplot(3,2,kk) 
    plot(x,v)
    line(xx,p)
    title(strcat('using h =', num2str(h)))
    E(kk) = norm(p - square(xx),2)*sqrt(h); %this is for demonstrating that the L^2 convergence is of order 1/2 
    error = abs(p - square(xx));
    error = error .* (error < 0.2); %the filter is arbitrary; we set anything with an error greater than 20% to zero.
    %one good reason for picking 20% is because we know that the error in
    %fourier series approximation with gibbs phenomenon is actually no more
    %than 16% 
    indices = find(error == max(error),1);
    locations(kk) = xx(indices); 
    overshoots(kk) = error(indices); 
end

figure %we need this because error_loglog produces a new graph and we want to see the old one 
error_loglog(dx,E)
error_table(dx,E)