% This program demonstrates a two-dimensional FDTD simulation(TM).
% The program mainly simulates the propagation of electromagnetic wave in 
% free space. The excitation source is a harmonic field source. The
% BPML boundary condition is used as the boundary condition.

% �ó�����ʾ�˶�άFDTDģ��(TM)��
% �ó�����Ҫģ���˵�Ų������ɿռ��еĴ���������ԴΪʱг��Դ���߽�����������
% BPML�߽�������

clc;
clear;
close all;

%***********************************************************************
% Fundamental constants ��������
%***********************************************************************

eps0 = 8.85e-12;	% permittivity of vacuum ��ս�糣��
mu0 = 4*pi*1e-7;	% permeability of vacuum ��մŵ���
c0 = 1/sqrt(mu0*eps0);	% speed of light ����
z0 = sqrt(mu0/eps0);    % Wave impedance of vacuum ����еĲ��迹

%***********************************************************************
% Mesh parameters �������
%***********************************************************************

Nx = 50;	% number of cells in 2D problem space ��ά����ռ��еĵ�Ԫ��
Ny = 50;
Nt = 300;	% number of iterations ��������
dx = 3e-2;	% space step �ռ䲽��
dy = 3e-2;
dt = 1/(2.0*c0*sqrt(1/dx^2+1/dy^2));	%time step ʱ�䲽��

%***********************************************************************
% Material properties ý������
%***********************************************************************

epsR = 1;	% relative permittivity ��Խ�糣��
muR = 1;	% relative permeability ��Դŵ���
sigE = 0;	% electric conductivity �絼��
sigH = 0;	% effective magnetism conductivity ��Ч�ŵ���

%***********************************************************************
% Boundary conditions �߽�����
%***********************************************************************

bpml = 8;
m = 4;
sigExmax = (m+1)/(sqrt(epsR)*150*pi*dx);
sigEymax = (m+1)/(sqrt(epsR)*150*pi*dy);
sigHxmax = (m+1)/(sqrt(epsR)*150*pi*dx);
sigHymax = (m+1)/(sqrt(epsR)*150*pi*dy);
sigEx = zeros(Nx+2*bpml+1,1);
sigEy = zeros(Ny+2*bpml+1,1);
sigHx = zeros(Nx+2*bpml,1);
sigHy = zeros(Ny+2*bpml,1);
for i = 1:bpml
    sigEx(bpml+1-i) = sigExmax*(i/bpml)^m; % ��
    sigEy(bpml+1-i) = sigEymax*(i/bpml)^m; % ��
    sigEx(Nx+bpml+1+i) = sigExmax*(i/bpml)^m;	% ��
    sigEy(Ny+bpml+1+i) = sigEymax*(i/bpml)^m;	% ��
    sigHx(bpml+1-i) = sigHxmax*(i/bpml)^m; % ��
    sigHy(bpml+1-i) = sigHymax*(i/bpml)^m; % ��
    sigHx(Nx+bpml+i) = sigHxmax*(i/bpml)^m;	% ��
    sigHy(Ny+bpml+i) = sigHymax*(i/bpml)^m;	% ��
end

%***********************************************************************
% Updating coefficients ����ϵ��
%***********************************************************************

CA = (eps0*epsR/dt-sigE/2)/(eps0*epsR/dt+sigE/2);
CB = 1/(eps0*epsR/dt+sigE/2);
CP = (mu0*muR/dt-sigH/2)/(mu0*muR/dt+sigH/2);
CQ = 1/(mu0*muR/dt+sigH/2);

%***********************************************************************
% Source excitation Դ����
%***********************************************************************

fre = 1.0e+9; % frequency ����ԴƵ��
Jx = round(Nx/2);	% position ����Դλ��
Jy = round(Ny/2);

%***********************************************************************
% Initializing field arrays ��ʼ����
%***********************************************************************

Hx = zeros(Nx+2*bpml+1,Ny+2*bpml);
Hy = zeros(Nx+2*bpml,Ny+2*bpml+1);
Ez = zeros(Nx+2*bpml+1,Ny+2*bpml+1);
Ezx = zeros(Nx+2*bpml+1,Ny+2*bpml+1);
Ezy = zeros(Nx+2*bpml+1,Ny+2*bpml+1);

%***********************************************************************
% BEGIN TIME-STEPPING LOOP ��ʼѭ��
%***********************************************************************

for n=1:Nt
    
    % Set excitation source ���ü���Դ
    Ez(Jx,Jy) = sin(2*pi*fre*n*dt);
    
    % Update magnetic field ���´ų�
    % Hx
    for i = 1+bpml:Nx+bpml+1
        for j = 1+bpml:Ny+bpml
            Hx(i,j) = CP*Hx(i,j)-CQ*(Ez(i,j+1)-Ez(i,j))/dy;
        end
    end
    % Hy
    for i = 1+bpml:Nx+bpml
        for j = 1+bpml:Ny+bpml+1
            Hy(i,j) = CP*Hy(i,j)+CQ*(Ez(i+1,j)-Ez(i,j))/dx;
        end
    end
    % PML
    % ��
    for i = 1:bpml
        for j = 1+bpml:Ny+bpml
            Hx(i,j) = Hx(i,j)-z0/2*(Ez(i,j+1)-Ez(i,j));
        end
    end
    for i = 1:bpml
        for j = 1+bpml:Ny+bpml+1
            Hy(i,j) = exp(-sigHx(i)*dt/mu0)*Hy(i,j)+(1-exp(-sigHx(i)*dt/mu0))/(dx*sigHx(i))*(Ez(i+1,j)-Ez(i,j));
        end
    end
    % ��
    for i = Nx+bpml+2:Nx+2*bpml+1
        for j = 1+bpml:Ny+bpml
            Hx(i,j) = Hx(i,j)-z0/2*(Ez(i,j+1)-Ez(i,j));
        end
    end
    for i = Nx+bpml+1:Nx+2*bpml
        for j = 1+bpml:Ny+bpml+1
            Hy(i,j) = exp(-sigHx(i)*dt/mu0)*Hy(i,j)+(1-exp(-sigHx(i)*dt/mu0))/(dx*sigHx(i))*(Ez(i+1,j)-Ez(i,j));
        end
    end
    % ��
    for i = 1+bpml:Nx+bpml+1
        for j = 1:bpml
            Hx(i,j) = exp(-sigHy(j)*dt/mu0)*Hx(i,j)-(1-exp(-sigHy(j)*dt/mu0))/(dy*sigHy(j))*(Ez(i,j+1)-Ez(i,j));
        end
    end
    for i = 1+bpml:Nx+bpml
        for j = 1:bpml
            Hy(i,j) = Hy(i,j)+z0/2*(Ez(i+1,j)-Ez(i,j));
        end
    end
    % ��
    for i = 1+bpml:Nx+bpml+1
        for j = Ny+bpml+1:Ny+2*bpml
            Hx(i,j) = exp(-sigHy(j)*dt/mu0)*Hx(i,j)-(1-exp(-sigHy(j)*dt/mu0))/(dy*sigHy(j))*(Ez(i,j+1)-Ez(i,j));
        end
    end
    for i = 1+bpml:Nx+bpml
        for j = Ny+bpml+2:Ny+2*bpml
            Hy(i,j) = Hy(i,j)+z0/2*(Ez(i+1,j)-Ez(i,j));
        end
    end
    % �ĸ���
    % ����
    for i = 1:bpml
        for j = 1:bpml
            Hx(i,j) = exp(-sigHy(j)*dt/mu0)*Hx(i,j)-(1-exp(-sigHy(j)*dt/mu0))/(dy*sigHy(j))*(Ez(i,j+1)-Ez(i,j));
        end
    end
    for i = 1:bpml
        for j = 1:bpml
            Hy(i,j) = exp(-sigHx(i)*dt/mu0)*Hy(i,j)+(1-exp(-sigHx(i)*dt/mu0))/(dx*sigHx(i))*(Ez(i+1,j)-Ez(i,j));
        end
    end
    % ����
    for i = 1:bpml
        for j = Ny+bpml+1:Ny+2*bpml
            Hx(i,j) = exp(-sigHy(j)*dt/mu0)*Hx(i,j)-(1-exp(-sigHy(j)*dt/mu0))/(dy*sigHy(j))*(Ez(i,j+1)-Ez(i,j));
        end
    end
    for i = 1:bpml
        for j = Ny+bpml+2:Ny+2*bpml+1
            Hy(i,j) = exp(-sigHx(i)*dt/mu0)*Hy(i,j)+(1-exp(-sigHx(i)*dt/mu0))/(dx*sigHx(i))*(Ez(i+1,j)-Ez(i,j));
        end
    end
    % ����
    for i = Nx+bpml+2:Nx+2*bpml+1
        for j = Ny+bpml+1:Ny+2*bpml
            Hx(i,j) = exp(-sigHy(j)*dt/mu0)*Hx(i,j)-(1-exp(-sigHy(j)*dt/mu0))/(dy*sigHy(j))*(Ez(i,j+1)-Ez(i,j));
        end
    end
    for i = Nx+bpml+1:Nx+2*bpml
        for j = Ny+bpml+2:Ny+2*bpml+1
            Hy(i,j) = exp(-sigHx(i)*dt/mu0)*Hy(i,j)+(1-exp(-sigHx(i)*dt/mu0))/(dx*sigHx(i))*(Ez(i+1,j)-Ez(i,j));
        end
    end
    % ����
    for i = Nx+bpml+2:Nx+2*bpml+1
        for j = 1:bpml
            Hx(i,j) = exp(-sigHy(j)*dt/mu0)*Hx(i,j)-(1-exp(-sigHy(j)*dt/mu0))/(dy*sigHy(j))*(Ez(i,j+1)-Ez(i,j));
        end
    end
    for i = Nx+bpml+1:Nx+2*bpml
        for j = 1:bpml
            Hy(i,j) = exp(-sigHx(i)*dt/mu0)*Hy(i,j)+(1-exp(-sigHx(i)*dt/mu0))/(dx*sigHx(i))*(Ez(i+1,j)-Ez(i,j));
        end
    end

%             Hx(i,j) = exp(-sigHy(j)*dt/mu0)*Hx(i,j)-(1-exp(-sigHy(j)*dt/mu0))/(dy*sigHy(j))*(Ez(i,j+1)-Ez(i,j));
%             Hy(i,j) = exp(-sigHx(i)*dt/mu0)*Hy(i,j)+(1-exp(-sigHx(i)*dt/mu0))/(dx*sigHx(i))*(Ez(i+1,j)-Ez(i,j));
%             Hx(i,j) = Hx(i,j)-z0/2*(Ez(i,j+1)-Ez(i,j));
%             Hy(i,j) = Hy(i,j)+z0/2*(Ez(i+1,j)-Ez(i,j));

    % Update electric field ���µ糡
    % Ez
    for i = 1+bpml:Nx+bpml+1
        for j = 1+bpml:Ny+bpml+1
           Ez(i,j) = CA*Ez(i,j)+CB* ...
                      ((Hy(i,j)-Hy(i-1,j))/dx-(Hx(i,j)-Hx(i,j-1))/dy);
        end
    end
    % PML
    % ��
    for i = 2:bpml
        for j = 1+bpml:Ny+bpml+1
            Ezx(i,j) = exp(-sigEx(i)*dt/eps0)*Ezx(i,j)+(1-exp(-sigEx(i)*dt/eps0))/(dx*sigEx(i))*(Hy(i,j)-Hy(i-1,j));
            Ezy(i,j) = Ezy(i,j)-1/(2*z0)*(Hx(i,j)-Hx(i,j-1));
            Ez(i,j) = Ezx(i,j)+Ezy(i,j);
        end
    end
    % ��
    for i = Nx+bpml+2:Nx+2*bpml
        for j = 1+bpml:Ny+bpml+1
            Ezx(i,j) = exp(-sigEx(i)*dt/eps0)*Ezx(i,j)+(1-exp(-sigEx(i)*dt/eps0))/(dx*sigEx(i))*(Hy(i,j)-Hy(i-1,j));
            Ezy(i,j) = Ezy(i,j)-1/(2*z0)*(Hx(i,j)-Hx(i,j-1));
            Ez(i,j) = Ezx(i,j)+Ezy(i,j);
        end
    end
    % ��
    for i = 1+bpml:Nx+bpml+1
        for j = 2:bpml
            Ezx(i,j) = Ezx(i,j)+1/(2*z0)*(Hy(i,j)-Hy(i-1,j));
            Ezy(i,j) = exp(-sigEy(j)*dt/eps0)*Ezy(i,j)-(1-exp(-sigEy(j)*dt/eps0))/(dy*sigEy(j))*(Hx(i,j)-Hx(i,j-1));
            Ez(i,j) = Ezx(i,j)+Ezy(i,j);
        end
    end
    % ��
    for i = 1+bpml:Nx+bpml+1
        for j = Ny+bpml+2:Ny+2*bpml
            Ezx(i,j) = Ezx(i,j)+1/(2*z0)*(Hy(i,j)-Hy(i-1,j));
            Ezy(i,j) = exp(-sigEy(j)*dt/eps0)*Ezy(i,j)-(1-exp(-sigEy(j)*dt/eps0))/(dy*sigEy(j))*(Hx(i,j)-Hx(i,j-1));
            Ez(i,j) = Ezx(i,j)+Ezy(i,j);
        end
    end
    % �ĸ���
    % ����
    for i = 2:bpml
        for j = 2:bpml
            Ezx(i,j) = exp(-sigEx(i)*dt/eps0)*Ezx(i,j)+(1-exp(-sigEx(i)*dt/eps0))/(dx*sigEx(i))*(Hy(i,j)-Hy(i-1,j));
            Ezy(i,j) = exp(-sigEy(j)*dt/eps0)*Ezy(i,j)-(1-exp(-sigEy(j)*dt/eps0))/(dy*sigEy(j))*(Hx(i,j)-Hx(i,j-1));
            Ez(i,j) = Ezx(i,j)+Ezy(i,j);
        end
    end
    % ����
    for i = 2:bpml
        for j = Ny+bpml+2:Ny+2*bpml
            Ezx(i,j) = exp(-sigEx(i)*dt/eps0)*Ezx(i,j)+(1-exp(-sigEx(i)*dt/eps0))/(dx*sigEx(i))*(Hy(i,j)-Hy(i-1,j));
            Ezy(i,j) = exp(-sigEy(j)*dt/eps0)*Ezy(i,j)-(1-exp(-sigEy(j)*dt/eps0))/(dy*sigEy(j))*(Hx(i,j)-Hx(i,j-1));
            Ez(i,j) = Ezx(i,j)+Ezy(i,j);
        end
    end
    % ����
    for i = Nx+bpml+2:Nx+2*bpml
        for j = Ny+bpml+2:Ny+2*bpml
            Ezx(i,j) = exp(-sigEx(i)*dt/eps0)*Ezx(i,j)+(1-exp(-sigEx(i)*dt/eps0))/(dx*sigEx(i))*(Hy(i,j)-Hy(i-1,j));
            Ezy(i,j) = exp(-sigEy(j)*dt/eps0)*Ezy(i,j)-(1-exp(-sigEy(j)*dt/eps0))/(dy*sigEy(j))*(Hx(i,j)-Hx(i,j-1));
            Ez(i,j) = Ezx(i,j)+Ezy(i,j);
        end
    end
    % ����
    for i = Nx+bpml+2:Nx+2*bpml
        for j = 2:bpml
            Ezx(i,j) = exp(-sigEx(i)*dt/eps0)*Ezx(i,j)+(1-exp(-sigEx(i)*dt/eps0))/(dx*sigEx(i))*(Hy(i,j)-Hy(i-1,j));
            Ezy(i,j) = exp(-sigEy(j)*dt/eps0)*Ezy(i,j)-(1-exp(-sigEy(j)*dt/eps0))/(dy*sigEy(j))*(Hx(i,j)-Hx(i,j-1));
            Ez(i,j) = Ezx(i,j)+Ezy(i,j);
        end
    end
%             Ezx(i,j) = exp(-sigEx(i)*dt/eps0)*Ezx(i,j)+(1-exp(-sigEx(i)*dt/eps0))/(dx*sigEx(i))*(Hy(i,j)-Hy(i-1,j));
%             Ezy(i,j) = Ezy(i,j)-1/(2*z0)*(Hx(i,j)-Hx(i,j-1));
%             Ez(i,j) = Ezx(i,j)+Ezy(i,j);
%             Ezx(i,j) = Ezx(i,j)+1/(2*z0)*(Hy(i,j)-Hy(i-1,j));
%             Ezy(i,j) = exp(-sigEy(j)*dt/eps0)*Ezy(i,j)-(1-exp(-sigEy(j)*dt/eps0))/(dy*sigEy(j))*(Hx(i,j)-Hx(i,j-1));
%             Ez(i,j) = Ezx(i,j)+Ezy(i,j);
         
    % Set boundary conditions ���ñ߽�����
    
    % Visualize fields ���ӻ���    
    imagesc(Ez');
    shading flat;caxis([-1.0 1.0]);axis image;axis xy; 
    title(['Ez, step ',int2str(n)]);xlabel('i');ylabel('j');
    drawnow;

end

%***********************************************************************
% END TIME-STEPPING LOOP ����ѭ��
%***********************************************************************