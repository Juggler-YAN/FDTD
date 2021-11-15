% This program demonstrates a two-dimensional FDTD simulation(Modified TE).
% The program mainly simulates the propagation of electromagnetic wave in 
% free space. The excitation source is a harmonic field source. The
% simplest PEC boundary condition is used as the boundary condition.

% �ó�����ʾ�˶�άFDTDģ��(�޸ĺ��TE)��
% �ó�����Ҫģ���˵�Ų������ɿռ��еĴ���������ԴΪʱг��Դ���߽�����������
% ��򵥵�PEC�߽�������

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

Nx = 50;	% number of cells in 2D problem space ��ά����ռ��еĵ�Ԫ��
Ny = 50;
Nt = 150;	% number of iterations ��������
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
Jx = round(Nx/2);	% position ����Դλ��
Jy = round(Ny/2);

%***********************************************************************
% Initializing field arrays ��ʼ����
%***********************************************************************

Ex = zeros(Nx,Ny+1);
Ey = zeros(Nx+1,Ny);
Hz = zeros(Nx+1,Ny+1);

%***********************************************************************
% BEGIN TIME-STEPPING LOOP ��ʼѭ��
%***********************************************************************

for n=1:Nt
    
    % Set excitation source ���ü���Դ
    Hz(Jx,Jy) = sin(2*pi*fre*n*dt);
    
    % Update electric field ���µ糡
    % Ex
    for i = 2:Nx
        for j = 1:Ny
            Ex(i,j) = CA*Ex(i,j)+CB*(Hz(i,j+1)-Hz(i,j))/dy;
        end
    end
    % Ey
    for i = 1:Nx
        for j = 2:Ny
            Ey(i,j) = CA*Ey(i,j)-CB*(Hz(i+1,j)-Hz(i,j))/dx;
        end
    end

    % Update magnetic field ���´ų�
    % Hz
    for i = 2:Nx
        for j = 2:Ny
           Hz(i,j) = CP*Hz(i,j)-CQ* ...
                      ((Ey(i,j)-Ey(i-1,j))/dx-(Ex(i,j)-Ex(i,j-1))/dy);
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