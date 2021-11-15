% This program demonstrates a two-dimensional FDTD simulation(Modified TE).
% The program mainly simulates the propagation of electromagnetic wave in 
% free space. The excitation source is a harmonic field source. The
% BPML boundary condition is used as the boundary condition.

% �ó�����ʾ�˶�άFDTDģ��(�޸ĺ��TE)��
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
sigEx = zeros(Nx+2*bpml,1);
sigEy = zeros(Ny+2*bpml,1);
sigHx = zeros(Nx+2*bpml+1,1);
sigHy = zeros(Ny+2*bpml+1,1);
for i = 1:bpml
    sigEx(bpml+1-i) = sigExmax*(i/bpml)^m; % ��
    sigEy(bpml+1-i) = sigEymax*(i/bpml)^m; % ��
    sigEx(Nx+bpml+i) = sigExmax*(i/bpml)^m;	% ��
    sigEy(Ny+bpml+i) = sigEymax*(i/bpml)^m;	% ��
    sigHx(bpml+1-i) = sigHxmax*(i/bpml)^m; % ��
    sigHy(bpml+1-i) = sigHymax*(i/bpml)^m; % ��
    sigHx(Nx+bpml+1+i) = sigHxmax*(i/bpml)^m;	% ��
    sigHy(Ny+bpml+1+i) = sigHymax*(i/bpml)^m;	% ��
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

Ex = zeros(Nx+2*bpml+1,Ny+2*bpml);
Ey = zeros(Nx+2*bpml,Ny+2*bpml+1);
Hz = zeros(Nx+2*bpml+1,Ny+2*bpml+1);
Hzx = zeros(Nx+2*bpml+1,Ny+2*bpml+1);
Hzy = zeros(Nx+2*bpml+1,Ny+2*bpml+1);

%***********************************************************************
% BEGIN TIME-STEPPING LOOP ��ʼѭ��
%***********************************************************************

for n=1:Nt
    
    % Set excitation source ���ü���Դ
    Hz(Jx,Jy) = sin(2*pi*fre*n*dt);
    
    % Update electric field ���µ糡
    % Ex
    for i = 1+bpml:Nx+bpml+1
        for j = 1+bpml:Ny+bpml
            Ex(i,j) = CA*Ex(i,j)+CB*(Hz(i,j+1)-Hz(i,j))/dy;
        end
    end
    % Ey
    for i = 1+bpml:Nx+bpml
        for j = 1+bpml:Ny+bpml+1
            Ey(i,j) = CA*Ey(i,j)-CB*(Hz(i+1,j)-Hz(i,j))/dx;
        end
    end
    % PML
    % ��
    for i = 1:bpml
        for j = 1+bpml:Ny+bpml
            Ex(i,j) = Ex(i,j)+z0/2*(Hz(i,j+1)-Hz(i,j));
        end
    end
    for i = 1:bpml
        for j = 1+bpml:Ny+bpml+1
            Ey(i,j) = exp(-sigEx(i)*dt/eps0)*Ey(i,j)-(1-exp(-sigEx(i)*dt/eps0))/(dx*sigEx(i))*(Hz(i+1,j)-Hz(i,j));
        end
    end
    % ��
    for i = Nx+bpml+2:Nx+2*bpml+1
        for j = 1+bpml:Ny+bpml
            Ex(i,j) = Ex(i,j)+z0/2*(Hz(i,j+1)-Hz(i,j));
        end
    end
    for i = Nx+bpml+1:Nx+2*bpml
        for j = 1+bpml:Ny+bpml+1
            Ey(i,j) = exp(-sigEx(i)*dt/eps0)*Ey(i,j)-(1-exp(-sigEx(i)*dt/eps0))/(dx*sigEx(i))*(Hz(i+1,j)-Hz(i,j));
        end
    end
    % ��
    for i = 1+bpml:Nx+bpml+1
        for j = 1:bpml
            Ex(i,j) = exp(-sigEy(j)*dt/eps0)*Ex(i,j)+(1-exp(-sigEy(j)*dt/eps0))/(dy*sigEy(j))*(Hz(i,j+1)-Hz(i,j));
        end
    end
    for i = 1+bpml:Nx+bpml
        for j = 1:bpml
            Ey(i,j) = Ey(i,j)-z0/2*(Hz(i+1,j)-Hz(i,j));
        end
    end
    % ��
    for i = 1+bpml:Nx+bpml+1
        for j = Ny+bpml+1:Ny+2*bpml
            Ex(i,j) = exp(-sigEy(j)*dt/eps0)*Ex(i,j)+(1-exp(-sigEy(j)*dt/eps0))/(dy*sigEy(j))*(Hz(i,j+1)-Hz(i,j));
        end
    end
    for i = 1+bpml:Nx+bpml
        for j = Ny+bpml+2:Ny+2*bpml+1
            Ey(i,j) = Ey(i,j)-z0/2*(Hz(i+1,j)-Hz(i,j));
        end
    end
    % �ĸ���
    % ����
    for i = 1:bpml
        for j = 1:bpml
            Ex(i,j) = exp(-sigEy(j)*dt/eps0)*Ex(i,j)+(1-exp(-sigEy(j)*dt/eps0))/(dy*sigEy(j))*(Hz(i,j+1)-Hz(i,j));
        end
    end
    for i = 1:bpml
        for j = 1:bpml
            Ey(i,j) = exp(-sigEx(i)*dt/eps0)*Ey(i,j)-(1-exp(-sigEx(i)*dt/eps0))/(dx*sigEx(i))*(Hz(i+1,j)-Hz(i,j));
        end
    end
    % ����
    for i = 1:bpml
        for j = Ny+bpml+1:Ny+2*bpml
            Ex(i,j) = exp(-sigEy(j)*dt/eps0)*Ex(i,j)+(1-exp(-sigEy(j)*dt/eps0))/(dy*sigEy(j))*(Hz(i,j+1)-Hz(i,j));
        end
    end
    for i = 1:bpml
        for j = Ny+bpml+2:Ny+2*bpml+1
            Ey(i,j) = exp(-sigEx(i)*dt/eps0)*Ey(i,j)-(1-exp(-sigEx(i)*dt/eps0))/(dx*sigEx(i))*(Hz(i+1,j)-Hz(i,j));
        end
    end
    % ����
    for i = Nx+bpml+2:Nx+2*bpml+1
        for j = Ny+bpml+1:Ny+2*bpml
            Ex(i,j) = exp(-sigEy(j)*dt/eps0)*Ex(i,j)+(1-exp(-sigEy(j)*dt/eps0))/(dy*sigEy(j))*(Hz(i,j+1)-Hz(i,j));
        end
    end
    for i = Nx+bpml+1:Nx+2*bpml
        for j = Ny+bpml+2:Ny+2*bpml+1
            Ey(i,j) = exp(-sigEx(i)*dt/eps0)*Ey(i,j)-(1-exp(-sigEx(i)*dt/eps0))/(dx*sigEx(i))*(Hz(i+1,j)-Hz(i,j));
        end
    end
    % ����
    for i = Nx+bpml+2:Nx+2*bpml+1
        for j = 1:bpml
            Ex(i,j) = exp(-sigEy(j)*dt/eps0)*Ex(i,j)+(1-exp(-sigEy(j)*dt/eps0))/(dy*sigEy(j))*(Hz(i,j+1)-Hz(i,j));
        end
    end
    for i = Nx+bpml+1:Nx+2*bpml
        for j = 1:bpml
            Ey(i,j) = exp(-sigEx(i)*dt/eps0)*Ey(i,j)-(1-exp(-sigEx(i)*dt/eps0))/(dx*sigEx(i))*(Hz(i+1,j)-Hz(i,j));
        end
    end

    % Update magnetic field ���´ų�
    % Hz
    for i = 1+bpml:Nx+bpml+1
        for j = 1+bpml:Ny+bpml+1
           Hz(i,j) = CP*Hz(i,j)-CQ* ...
                      ((Ey(i,j)-Ey(i-1,j))/dx-(Ex(i,j)-Ex(i,j-1))/dy);
        end
    end
    % PML
    % ��
    for i = 2:bpml
        for j = 1+bpml:Ny+bpml+1
            Hzx(i,j) = exp(-sigHx(i)*dt/mu0)*Hzx(i,j)-(1-exp(-sigHx(i)*dt/mu0))/(dx*sigHx(i))*(Ey(i,j)-Ey(i-1,j));
            Hzy(i,j) = Hzy(i,j)+1/(2*z0)*(Ex(i,j)-Ex(i,j-1));
            Hz(i,j) = Hzx(i,j)+Hzy(i,j);
        end
    end
    % ��
    for i = Nx+bpml+2:Nx+2*bpml
        for j = 1+bpml:Ny+bpml+1
            Hzx(i,j) = exp(-sigHx(i)*dt/mu0)*Hzx(i,j)-(1-exp(-sigHx(i)*dt/mu0))/(dx*sigHx(i))*(Ey(i,j)-Ey(i-1,j));
            Hzy(i,j) = Hzy(i,j)+1/(2*z0)*(Ex(i,j)-Ex(i,j-1));
            Hz(i,j) = Hzx(i,j)+Hzy(i,j);
        end
    end
    % ��
    for i = 1+bpml:Nx+bpml+1
        for j = 2:bpml
            Hzx(i,j) = Hzx(i,j)-1/(2*z0)*(Ey(i,j)-Ey(i-1,j));
            Hzy(i,j) = exp(-sigHy(j)*dt/mu0)*Hzy(i,j)+(1-exp(-sigHy(j)*dt/mu0))/(dy*sigHy(j))*(Ex(i,j)-Ex(i,j-1));
            Hz(i,j) = Hzx(i,j)+Hzy(i,j);
        end
    end
    % ��
    for i = 1+bpml:Nx+bpml+1
        for j = Ny+bpml+2:Ny+2*bpml
            Hzx(i,j) = Hzx(i,j)-1/(2*z0)*(Ey(i,j)-Ey(i-1,j));
            Hzy(i,j) = exp(-sigHy(j)*dt/mu0)*Hzy(i,j)+(1-exp(-sigHy(j)*dt/mu0))/(dy*sigHy(j))*(Ex(i,j)-Ex(i,j-1));
            Hz(i,j) = Hzx(i,j)+Hzy(i,j);
        end
    end
    % �ĸ���
    % ����
    for i = 2:bpml
        for j = 2:bpml
            Hzx(i,j) = exp(-sigHx(i)*dt/mu0)*Hzx(i,j)-(1-exp(-sigHx(i)*dt/mu0))/(dx*sigHx(i))*(Ey(i,j)-Ey(i-1,j));
            Hzy(i,j) = exp(-sigHy(j)*dt/mu0)*Hzy(i,j)+(1-exp(-sigHy(j)*dt/mu0))/(dy*sigHy(j))*(Ex(i,j)-Ex(i,j-1));
            Hz(i,j) = Hzx(i,j)+Hzy(i,j);
        end
    end
    % ����
    for i = 2:bpml
        for j = Ny+bpml+2:Ny+2*bpml
            Hzx(i,j) = exp(-sigHx(i)*dt/mu0)*Hzx(i,j)-(1-exp(-sigHx(i)*dt/mu0))/(dx*sigHx(i))*(Ey(i,j)-Ey(i-1,j));
            Hzy(i,j) = exp(-sigHy(j)*dt/mu0)*Hzy(i,j)+(1-exp(-sigHy(j)*dt/mu0))/(dy*sigHy(j))*(Ex(i,j)-Ex(i,j-1));
            Hz(i,j) = Hzx(i,j)+Hzy(i,j);
        end
    end
    % ����
    for i = Nx+bpml+2:Nx+2*bpml
        for j = Ny+bpml+2:Ny+2*bpml
            Hzx(i,j) = exp(-sigHx(i)*dt/mu0)*Hzx(i,j)-(1-exp(-sigHx(i)*dt/mu0))/(dx*sigHx(i))*(Ey(i,j)-Ey(i-1,j));
            Hzy(i,j) = exp(-sigHy(j)*dt/mu0)*Hzy(i,j)+(1-exp(-sigHy(j)*dt/mu0))/(dy*sigHy(j))*(Ex(i,j)-Ex(i,j-1));
            Hz(i,j) = Hzx(i,j)+Hzy(i,j);
        end
    end
    % ����
    for i = Nx+bpml+2:Nx+2*bpml
        for j = 2:bpml
            Hzx(i,j) = exp(-sigHx(i)*dt/mu0)*Hzx(i,j)-(1-exp(-sigHx(i)*dt/mu0))/(dx*sigHx(i))*(Ey(i,j)-Ey(i-1,j));
            Hzy(i,j) = exp(-sigHy(j)*dt/mu0)*Hzy(i,j)+(1-exp(-sigHy(j)*dt/mu0))/(dy*sigHy(j))*(Ex(i,j)-Ex(i,j-1));
            Hz(i,j) = Hzx(i,j)+Hzy(i,j);
        end
    end
         
    % Set boundary conditions ���ñ߽�����
    
    % Visualize fields ���ӻ���    
    imagesc(Hz');
    shading flat;caxis([-1.0 1.0]);axis image;axis xy; 
    title(['Hz, step ',int2str(n)]);xlabel('i');ylabel('j');
    drawnow;

end

%***********************************************************************
% END TIME-STEPPING LOOP ����ѭ��
%***********************************************************************