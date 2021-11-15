% This program demonstrates a one-dimensional FDTD simulation.
% The program mainly simulates the propagation of electromagnetic wave in Z 
% directionin free space. The excitation source is a harmonic field source.
% The Mur absorbing boundary condition is used as the boundary condition.

% �ó�����ʾ��һάFDTDģ�⡣
% �ó�����Ҫģ���˵�Ų������ɿռ���z�����ϵĴ���������ԴΪʱг��Դ���߽�����
% ������Mur���ձ߽�������

clc;
clear;
close all;

%***********************************************************************
% Fundamental constants ��������
%***********************************************************************

eps0 = 8.85e-12;	% permittivity of vacuum ��ս�糣��
mu0 = 4*pi*1e-7;	% permeability of vacuum ��մŵ���
c0 = 1/sqrt(mu0*eps0);	% speed of light ����

%***********************************************************************
% Mesh parameters �������
%***********************************************************************

Nz = 100;	% number of cells in 1D problem space һά����ռ��еĵ�Ԫ��
Nt = 150;	% number of iterations ��������
dz = 3e-2;	% space step �ռ䲽��
dt = 1/(2.0*c0*sqrt(1/dz^2));	%time step ʱ�䲽��

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

% PEC�߽�����

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
Jz = round(Nz/2);   % position ����Դλ��

%***********************************************************************
% Initializing field arrays ��ʼ����
%***********************************************************************

Ex = zeros(Nz+1,1);
Hy = zeros(Nz,1);

%***********************************************************************
% BEGIN TIME-STEPPING LOOP ��ʼѭ��
%***********************************************************************

for n=1:Nt
    
    % Set excitation source ���ü���Դ
%     Ex(1) = sin(2.0*pi*fre*dt*n);
%     Ex(Jz) = 1;
    Ex(Jz) = sin(2.0*pi*fre*dt*n);
%     Ex(Nz+1) = sin(2.0*pi*fre*dt*n);

    % Update magnetic field ���´ų�
    for k = 1:Nz
        Hy(k) = CP*Hy(k)-CQ*(Ex(k+1)-Ex(k))/dz;
    end
    
    tempx = Ex;
    
    % Update electric field ���µ糡
    for k = 2:Nz
        Ex(k) = CA*Ex(k)-CB*(Hy(k)-Hy(k-1))/dz;
    end
   
    % Set boundary conditions ���ñ߽�����
    % Murһ�׽��Ʒ�
%     Ex(1) = tempx(2)+((c0*dt-dz)/(c0*dt+dz))*(Ex(2)-tempx(1));
%     Ex(Nz+1) = tempx(Nz)+((c0*dt-dz)/(c0*dt+dz))*(Ex(Nz)-tempx(Nz+1));
    % �в���ʱ��
%     Ex(1) = tempx(1)-c0*dt/dz*(tempx(1)-tempx(2));
%     Ex(Nz+1) = tempx(Nz+1)-c0*dt/dz*(tempx(Nz+1)-tempx(Nz));
    % ���迹��
    Z0 = sqrt((muR*mu0)/(epsR*eps0));
    Ex(1) = -Z0*Hy(1);
    Ex(Nz+1) = Z0*Hy(Nz);
    
    % Visualize fields ���ӻ���
    subplot(2,1,1),plot(Ex,'r');
    xlabel('z');ylabel(['Ex, step ',int2str(n)]);axis([0 100 -2 2]);
    subplot(2,1,2),plot(Hy,'b');
    xlabel('z');ylabel(['Hy, step ',int2str(n)]);axis([0 100 -5e-3 5e-3]);
    drawnow;

end

%***********************************************************************
% END TIME-STEPPING LOOP ����ѭ��
%***********************************************************************