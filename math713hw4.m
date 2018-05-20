%% hw 3 
%% problem 12.7
close all
clear 
clc

%demonstration of code working
N = 10000; %N = 10^5
f = @(z) log(1 + 0.5*z); 

jj = 1:N; %the index of the taylor coefficient
theta = jj*2*pi/N; 
g = f(exp(1i*theta)); 
a_n = 1/N*fft(g); 
first_twenty = real(a_n(1:20)); 

for ii = 1:20
    fprintf('The %dth coefficient of the Taylor series expansion is %0.5f\n',ii-1,first_twenty(ii))
end

%asymptotic convergence 
L = 15; %number of trials we will run 
N_eval = logspace(2,5,L);
E = zeros(1,L-1); 
previous = zeros(20,1);
count = 1; 

for N = N_eval
    jj = 1:N; %the index of the taylor coefficient
    theta = jj*2*pi/N; 
    g = f(exp(1i*theta)); 
    a_n = 1/N*fft(g); 
    first_twenty = real(a_n(1:20)); 
    error = norm(first_twenty - previous,inf); %compute the error 
    E(count) = error; %store the error 
    
    %update 
    previous = first_twenty; 
    count = count + 1; 
end

loglog(N_eval,E)
title('plot of error vs. points used')
xlabel('N, number of points used')
ylabel('Error')

%% 11.1 
%eigenvalues of Laplacian on annulus from analytical solution 
close all
clear 
clc 

eqn = @(alpha) besselj(0,alpha).*bessely(0,2*alpha) - bessely(0,alpha).*besselj(0,2*alpha);
max = 40;
xx = -max:0.001:max;
yy = eqn(xx); 
yyzero = 0*xx; 
plot(xx,yy,'b',xx,yyzero,'r')
title('J_0(\nu)Y_0(2\nu) - J_0(2\nu)Y_0(\nu) = 0')
xlabel('\nu')

test = 1:60; 
nu_crit = fsolve(eqn,test); 

alpha = zeros(11,1); 
alpha(1) = nu_crit(1); 
do_nothing = 0;

count = 1; 
for ii = 1:60
    %determination of whether or not the solutions in holla are close
    if abs(nu_crit(ii) - nu_crit(ii+1)) < 1e-6 %i.e. the two consecutive entries are sufficiently close
        do_nothing = do_nothing + 1;  
    else %they are not sufficiently close 
        alpha(count + 1) = nu_crit(ii + 1); 
        count = count + 1; 
    end
    
    %termination of search 
    if count > 10
        break
    end
end

hold on 
plot(alpha,zeros(size(alpha)),'kx')
alpha = alpha.^2; %fam contains the eigenvalues of the laplacian on the annulus
alpha = alpha(1:10); %to make sure we only get the 10

%eigenvalues of Laplacian from Chebyshev matrices
N = 28; 
[D,x] = cheb(N); 
r = x/2 + 3/2; %define r on [1,2]
D = 2*D; %need to rescale the differentiation matrix. 
D2 = D^2; 
D1 = D;
D1r = diag(1./r) * D1; 

D2 = D2(2:N,2:N); 
D1r = D1r(2:N,2:N); 

L = D2 + D1r; 
L = -L; 
eigen = sort(eig(L));
eigen10 = eigen(1:10); 

for ii = 1:10
    fprintf('The eigenvalue as computed from solving the Sturm-Liouville problem is %4.8f \n The eigenvalue as computed by Chebyshev matrices is %4.8f \n',alpha(ii),eigen10(ii))
end

%% problem 10.4
close all
clear all 
clc 

N = 26; 
[D,x] = cheb(N); 
D2 = D^2; 
D2 = D2(2:N,2:N); 

F = @(u) D2*u + exp(u); 
dx = abs(x(2:end) - x(1:end-1)); 
hmin = min(dx); 
dt = 1/10*hmin^2; 
zeroindex = find(abs(x) < 1e-6); 

u = zeros(N-1,1); %u is a column of zeros; 
M = ceil(3.5/dt); %maximum number of iterations
for ii = 1:3*M
    k1 = dt*F(u);
    k2 = dt*F(u + k1/2); 
    k3 = dt*F(u + k2/2);
    k4 = dt*F(u + k3); 
    u = u + 1/6*(k1 + 2*k2 + 2*k3 + k4);
    u0 = u(zeroindex);
    if ii == M
        heat = u0; 
    elseif abs(u0 - 5) < 1e-3
        time = ii*dt;
        break
    end
    
end
 
parta_answer = heat; 

fprintf('u(0,3.5) = %2.8f \n The time t for which u(0,t) = 5 is %2.8f\n',parta_answer,time) 
