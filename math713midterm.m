%% midterm 7.4
close all
clc
clear all

values = 6:2:60;
finally = zeros(size(values));

for iii = 1:length(values)
  
hello = cputime;
  N = values(iii);
  [D,x] = cheb(N); D2 = D^2; D2 = D2(2:N,2:N);
  u = zeros(N-1,1);
  change = 1; it = 0;
  while change > 1e-15                   % fixed-point iteration 
%     unew = D2\exp(u); 
%     change = norm(unew-u,inf);
%     u = unew; it = it+1;
    deriv = D2 - diag(exp(u)); 
    increm = deriv\(D2*u - exp(u)); 
    unew = u - increm; 
    change = norm(unew - u, inf); 
    u = unew; 
    it = it + 1; 
  end
  u = [0;u;0];
  clf, subplot('position',[.1 .4 .8 .5])
  plot(x,u,'.','markersize',16)
  xx = -1:.01:1;
  uu = polyval(polyfit(x,u,N),xx);
  line(xx,uu), grid on
  title(sprintf('no. steps = %d      u(0) =%18.14f',it,u(N/2+1)))
 goodbye = cputime; 
 format long
 interval_of_time(iii) = goodbye - hello;
 finally(iii) = u(N/2+1);
end
hh = 1./values; 
roc1 = (finally(2:end) - u(N/2+1))./(finally(1:end-1) - u(N/2+1)).^2;
figure
stem(roc1)
title('estimate of rate of convergence')
figure
plot(values,interval_of_time,'-o')
xlabel('N')
ylabel('interval of time')
title('N vs time')

for jjj = 1:length(values)
    disp(sprintf('N = %d, u(0) = %2.15f',values(jjj),finally(jjj)))
end
  %% problem 8.4
  
close all
clear all
clc

hello = cputime;
N = 24; 
[D,x] = cheb(N); 
y = x'; 
  dt = 6/N^2;
  D2 = D^2; 
  D2 = D2(2:N,2:N); %since we will only need the interior points. note: this is a matrix that is N - 1 x N - 1
  I = eye(N-1);
  L = kron(D2,I) + kron(I,D2);
  [xx,yy] = meshgrid(x,y); % xx, yy help establish initial condition 
  plotgap = round((1/3)/dt); dt = (1/3)/plotgap;
  vv = exp(-40*((xx-.4).^2 + yy.^2)); %this is the inital condition 
  vvint = vv(2:N,2:N);
  current = reshape(vvint,(N-1)^2,1);
  vvold = vv; %because if du/dt = 0, this means u(x,y,t = 0) = u(x,y,t = dt)
  old = current; 

% Time-stepping by leap frog formula:
  [ay,ax] = meshgrid([.56 .06],[.1 .55]); clf
  for n = 0:3*plotgap
    t = n*dt;
    if rem(n+.5,plotgap)<1     % plots at multiples of t=1/3
      i = n/plotgap+1;
      subplot('position',[ax(i) ay(i) .36 .36])
      [xxx,yyy] = meshgrid(-1:1/16:1,-1:1/16:1); %the plotting grid
      vvv = interp2(xx,yy,vv,xxx,yyy,'cubic');
      mesh(xxx,yyy,vvv), axis([-1 1 -1 1 -0.15 1])
      colormap(1e-6*[1 1 1]); title(['t = ' num2str(t)]), drawnow
    end
    advanced = (-old/dt^2 + (L + 2/dt^2* eye((N-1)^2) )* current) * dt^2; 
    old = current; %what was once current is now old
    current = advanced; %what was once advanced is now current 
    vvint  = reshape(current,N-1,N-1);
    vv(2:N,2:N) = vvint;
  
  end
  
  goodbye = cputime; 
  interval_of_time = goodbye - hello;
  
fprintf('the amount of time it takes for midterm project to run is %4.15f seconds\n',interval_of_time)
  %% problem 8.5
  
  close all
clear all
clc

hello = cputime;

N = 24; 
[D,x] = cheb(N); 
y = x'; 
  dt = 6/N^2;
  D2 = D^2; 
  D2 = D2(2:N,2:N); %since we will only need the interior points. note: this is a matrix that is N - 1 x N - 1
  I = eye(N-1);
  L = kron(D2,I) + kron(I,D2);
  [xx,yy] = meshgrid(x,y); % xx, yy help establish initial condition 
  plotgap = round((1/3)/dt); dt = (1/3)/plotgap;
  vv = exp(-40*((xx-.4).^2 + yy.^2)); %this is the inital condition 
  vvint = vv(2:N,2:N);
  current = reshape(vvint,(N-1)^2,1);
  vvold = vv; %because if du/dt = 0, this means u(x,y,t = 0) = u(x,y,t = dt)
  old = current;
  u0 = old; %this is the initial condition
  
    s0 = zeros(2*(N-1)^2); %this is the vector of initial conditions
    s0(1:(N-1)^2) = u0; %the inital condition on u
    
    TT = zeros(2*(N-1)^2);
    TT((N-1)^2+1:end,1:(N-1)^2) = L;  
    TT(1:(N-1)^2,(N-1)^2+1:end) = eye((N-1)^2);

% Time-stepping by leap frog formula:
  [ay,ax] = meshgrid([.56 .06],[.1 .55]); clf
  for n = 0:3*plotgap
    t = n*dt;
    
    if rem(n+.5,plotgap)<1     % plots at multiples of t=1/3
      i = n/plotgap+1;
      subplot('position',[ax(i) ay(i) .36 .36])
      [xxx,yyy] = meshgrid(-1:1/16:1,-1:1/16:1); %the plotting grid
      vvv = interp2(xx,yy,vv,xxx,yyy,'cubic');
      mesh(xxx,yyy,vvv), axis([-1 1 -1 1 -0.15 1])
      colormap(1e-6*[1 1 1]); title(['t = ' num2str(t)]), drawnow
    end
%     advanced = (-old/dt^2 + (L + 2/dt^2* eye((N-1)^2) )* current) * dt^2; 
%     old = current; %what was once current is now old
%     current = advanced; %what was once advanced is now current 
%     vvint  = reshape(current,N-1,N-1);
%     vv(2:N,2:N) = vvint;
%   then using expm
    % in essence TT = (0 I; A 0) 
    deliverance = expm(TT*t)*s0; 
    vvint = reshape(deliverance(1:(N-1)^2), N-1, N-1);
    vv(2:N,2:N) = vvint;
  end

  goodbye = cputime;
  interval_of_time = goodbye - hello;  
fprintf('the amount of time it takes for midterm project matrix exponential to run is %4.15f seconds\n',interval_of_time)