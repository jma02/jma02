function dno_series

% Computes the DNO in series form by the FFT

clc; clear all;

%fid = fopen('dno.out','w');

% numerical parameters
%%%%%%%%%%%%%%%%%%%%%%

%Nx = input('Enter number of grid points Nx ==> ');
%Nx = 256;
Nx = 16384;

%L = input('Enter length of domain L ==> ');
%L = 2*pi;
L = 334*2*pi;

dx = L/Nx;
dk = 2*pi/L;
x = dx*[0:Nx-1]';
k = dk*[0:Nx/2,1-Nx/2:-1]';

% physical parameters
%%%%%%%%%%%%%%%%%%%%%

g = 1;

%h = input('Enter water depth h ==> ');
h = 1;

n0 = input('Enter wavenumber n0 (integer) ==> ');
%n0 = 1;

k0 = n0*dk;

%M = input('Enter truncation order M (integer) ==> ');   
M = 12;

mter = [0:M]';

% prescribed data
%%%%%%%%%%%%%%%%%

a0 = input('Enter surface amplitude a0 (< 1) ==> ');  
%a0 = 0.1;

om = sqrt(g*k0*tanh(k0*h));

eta = a0*cos(k0*x);
xi = a0*g*cosh(k0*(eta + h)).*sin(k0*x)./(om*cosh(k0*h));

eta_x = myifft(1i*k.*myfft(eta,Nx));

aux = k0*sinh(k0*(eta + h)).*sin(k0*x) - k0*eta_x.*cosh(k0*(eta + h)).*cos(k0*x);
Gexac = a0*g*aux./(om*cosh(k0*h));

Gerr = zeros(M+1,1);

G0 = k.*tanh(h*k);

for j=1:M+1
    fprintf(['m = %d',' terms\n'],j-1);
    tic
    Gnume = dno(eta,xi,k,G0,Nx,j-1);
    toc
    Gerr(j) = norm(Gnume - Gexac)/norm(Gexac);
end

semilogy(mter,Gerr,'-o');
xlabel('$M$','interpreter','latex','FontSize',16);
ylabel('Error','interpreter','latex','FontSize',16);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [G] = dno(eta,xi,k,G0,Nx,M)
% Sums up all terms to evaluate the DNO series

etam = zeros(Nx,M+1);

etam(:,0+1) = ones(Nx,1);
for m=1:M
    etam(:,m+1) = multiply(eta,etam(:,m-1+1),Nx,M)./m;
end

fxi = myfft(xi,Nx);
xi_x = myifft(1i*k.*fxi);

Gm = zeros(Nx,M+1); 
[Gm] = dno_n(etam,xi_x,fxi,k,G0,Nx,M);

G = Gm(:,0+1);
for m=1:M
    G = G + Gm(:,m+1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Gm] = dno_n(etam,xi_x,fxi,k,G0,Nx,M)
% Computes each individual term in the DNO series

Gm(:,0+1) = myifft(G0.*fxi);

for m=1:M
    if (mod(m,2) == 0)
        r = round(m/2);
        tmp = multiply(etam(:,m+1),xi_x,Nx,M);
        Gm(:,m+1) = -myifft(G0.*k.^(2*(r-1)).*(1i*k).*myfft(tmp,Nx));
        for s=0:r-1
            tmp = multiply(etam(:,2*(r-s)+1),Gm(:,2*s+1),Nx,M);
            Gm(:,m+1) = Gm(:,m+1) - myifft(k.^(2*(r-s)).*myfft(tmp,Nx));
            tmp = multiply(etam(:,(2*(r-s)-1)+1),Gm(:,(2*s+1)+1),Nx,M);
            Gm(:,m+1) = Gm(:,m+1) - myifft(G0.*k.^(2*(r-s-1)).*myfft(tmp,Nx));
        end
    else
        r = round((m+1)/2);
        tmp = multiply(etam(:,m+1),xi_x,Nx,M);
        Gm(:,m+1) = -myifft(k.^(2*(r-1)).*(1i*k).*myfft(tmp,Nx));
        for s=0:r-2
            tmp = multiply(etam(:,(2*(r-s)-1)+1),Gm(:,2*s+1),Nx,M);
            Gm(:,m+1) = Gm(:,m+1) - myifft(G0.*k.^(2*(r-s-1)).*myfft(tmp,Nx));
            tmp = multiply(etam(:,2*(r-s-1)+1),Gm(:,(2*s+1)+1),Nx,M);
            Gm(:,m+1) = Gm(:,m+1) - myifft(k.^(2*(r-s-1)).*myfft(tmp,Nx));
        end
        tmp = multiply(etam(:,1+1),Gm(:,(2*(r-1))+1),Nx,M);
        Gm(:,m+1) = Gm(:,m+1) - myifft(G0.*myfft(tmp,Nx));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [y] = multiply(a,b,Nx,M)
% Multiplies in physical space after Fourier extension

fac = 8;   % 1 or power of 2

Ny = fac*Nx;   % enlarged size for zero-padding in Fourier space

fa = myfft(a,Nx);
fb = myfft(b,Nx);

fya = zeros(Ny,1);
fyb = zeros(Ny,1);

fya(1:Nx/2+1) = fa(1:Nx/2+1);
fyb(1:Nx/2+1) = fb(1:Nx/2+1);

fya(Ny+1-Nx/2:Ny) = fa(Nx/2+1:Nx);
fyb(Ny+1-Nx/2:Ny) = fb(Nx/2+1:Nx);

ya = myifft(fya);
yb = myifft(fyb);

fw = myfft(ya.*yb,Nx);

fy = zeros(Nx,1);

fy(1:Nx/2+1) = fw(1:Nx/2+1);
fy(Nx/2+1:Nx) = fw(Ny+1-Nx/2:Ny);

y = fac*myifft(fy);

%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fy] = myfft(y,Nx)
% Direct Fourier transform

fy = fft(y);
fy(Nx/2+1) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%

function [y] = myifft(fy)
% Inverse Fourier transform for real functions

y = real(ifft(fy));
