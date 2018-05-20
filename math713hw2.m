%% math 713 hw 
%% problem 4.7
clear all
close all
clc

nn = zeros(1,16); %this will hold the information on number of gridpoints
eigen20 = zeros(20,16); %this will hold the actual numerical value of the eigenvalues

format long, format compact
  L = 8;                             % domain is [-L L], periodic
  count = 1; 
  for N = 20:2:100
    h = 2*pi/N; x = h*(1:N); x = L*(x-pi)/pi;
    column = [-pi^2/(3*h^2)-1/6 ...
       -.5*(-1).^(1:N-1)./sin(h*(1:N-1)/2).^2];
    D2 = (pi/L)^2*toeplitz(column);  % 2nd-order differentiation
    eigenvalues = sort(eig(-D2 + diag(x.^4)));
    N, eigenvalues
    nn(count) = N; 
    eigen20(:,count) = eigenvalues(1:20);
    count = count + 1; 
  end

  subplot(2,2,1)
  plot(nn,eigen20)
  title('values of the eigenvalues')
  xlabel('number of gridpoints')
  ylabel('eigenvalue')
  
  subplot(2,2,2)
  plot(nn(end-10:end),eigen20(:,end-10:end))
  title('values of the eigenvalues over HIGH N')
  xlabel('number of gridpoints')
  ylabel('eigenvalue')
  
  subplot(2,2,3)
  plot(nn,eigen20(1,:))
  title('value of the first eigenvalue')
  xlabel('number of gridpoints')
  ylabel('eigenvalue')
  
  subplot(2,2,4)
  plot(nn,eigen20(end,:))
  title('value of the last eigenvalue')
  xlabel('number of gridpoints')
  ylabel('eigenvalue')
  
  %% problem 5.1
  clear all
  close all
  clc
  
  n = linspace(4,18,6);
  n = floor(n); %to make sure that we get an integer
  for ii = 1:length(n)
      N = n(ii);
    xx = -1.01:.005:1.01; clf
  for i = 1:2
    if i==1, s = 'equispaced points'; x = -1 + 2*(0:N)/N; end
    if i==2, s = 'Chebyshev points';  x = cos(pi*(0:N)/N); end
%     subplot(2,2,i)
    u = 1./(1+16*x.^2); %this is the value of u at the gridpoints x_j 
    uu = 1./(1+16*xx.^2); %this is for plotting
    p = polyfit(x,u,N);              % interpolation (simple polynomial interpolation)
    pp = polyval(p,xx);              % evaluation of interpolant
%     plot(x,u,'.','markersize',13)
%     line(xx,pp)
%     axis([-1.1 1.1 -1 1.5]), title(s)
    error(ii,i) = norm(uu-pp,inf); %the column indicates whether or not equispaced or chebyshev was used; row indicates N 
%     text(-.5,-.5,['max error = ' num2str(error)])
  end
  end
  
  subplot(2,1,1)
  semilogy(n,error(:,1))
  title('error using equispaced points')
  
  subplot(2,1,2)
  semilogy(n,error(:,2))
  title('error using chebyshev points')