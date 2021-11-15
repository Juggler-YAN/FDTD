% This program demonstrates a two-dimensional FDTD simulation(Modified TE).
% The program mainly simulates the propagation of electromagnetic wave in 
% free space. The excitation source is a harmonic field source. The Mur
% absorbing boundary condition is used as the boundary condition.

% �ó�����ʾ�˶�άFDTDģ��(�޸ĺ��TE)��
% �ó�����Ҫģ���˵�Ų������ɿռ��еĴ���������ԴΪʱг��Դ���߽�����������
% Mur���ձ߽�������

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

Ex = zeros(Nx+1,Ny);
Ey = zeros(Nx,Ny+1);
Hz = zeros(Nx+1,Ny+1);

%***********************************************************************
% BEGIN TIME-STEPPING LOOP ��ʼѭ��
%***********************************************************************

temp2z = Hz;
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

    tempz = Hz;
    
    % Update magnetic field ���´ų�
    % Hz
    for i = 2:Nx
        for j = 2:Ny
           Hz(i,j) = CP*Hz(i,j)-CQ* ...
                      ((Ey(i,j)-Ey(i-1,j))/dx-(Ex(i,j)-Ex(i,j-1))/dy);
        end
    end 
    
    % Set boundary conditions ���ñ߽�����
    % ���
    % һ��Mur
%     % ��߽�
%     for j = 2:Ny
%         Hz(1,j) = tempz(2,j)+((c0*dt-dx)/(c0*dt+dx))*(Hz(2,j)-tempz(1,j));
%     end
%     % �ұ߽�
%     for j = 2:Ny
%         Hz(Nx+1,j) = tempz(Nx,j)+((c0*dt-dx)/(c0*dt+dx))*(Hz(Nx,j)-tempz(Nx+1,j));
%     end
%     % �±߽�
%     for i = 2:Nx
%         Hz(i,1) = tempz(i,2)+((c0*dt-dy)/(c0*dt+dy))*(Hz(i,2)-tempz(i,1));
%     end
%     % �ϱ߽�
%     for i = 2:Nx
%         Hz(i,Ny+1) = tempz(i,Ny)+((c0*dt-dy)/(c0*dt+dy))*(Hz(i,Ny)-tempz(i,Ny+1));
%     end
    % ����Mur
    % ��߽�
    for j = 2:Ny
        Hz(1,j) = tempz(2,j)+((c0*dt-dx)/(c0*dt+dx))*(Hz(2,j)-tempz(1,j))+ ...
                  (c0*c0*eps0*epsR*dt)/(2*(c0*dt+dx))*(dx/dy)*(Ex(1,j)-Ex(1,j-1)+Ex(2,j)-Ex(2,j-1));
    end
    % �ұ߽�
    for j = 2:Ny
        Hz(Nx+1,j) = tempz(Nx,j)+((c0*dt-dx)/(c0*dt+dx))*(Hz(Nx,j)-tempz(Nx+1,j))+ ...
                     (c0*c0*eps0*epsR*dt)/(2*(c0*dt+dx))*(dx/dy)*(Ex(Nx+1,j)-Ex(Nx+1,j-1)+Ex(Nx,j)-Ex(Nx,j-1));
    end
    % �±߽�
    for i = 2:Nx
        Hz(i,1) = tempz(i,2)+((c0*dt-dy)/(c0*dt+dy))*(Hz(i,2)-tempz(i,1))- ...
                  (c0*c0*eps0*epsR*dt)/(2*(c0*dt+dy))*(dy/dx)*(Ey(i,1)-Ey(i-1,1)+Ey(i,2)-Ey(i-1,2));
    end
    % �ϱ߽�
    for i = 2:Nx
        Hz(i,Ny+1) = tempz(i,Ny)+((c0*dt-dy)/(c0*dt+dy))*(Hz(i,Ny)-tempz(i,Ny+1))- ...
                     (c0*c0*eps0*epsR*dt)/(2*(c0*dt+dy))*(dy/dx)*(Ey(i,Ny+1)-Ey(i-1,Ny+1)+Ey(i,Ny)-Ey(i-1,Ny));
    end
    % �ǵ㣨����һ��
%     % ���½�
%     Hz(1,1)=tempz(2,2)+((c0*dt-sqrt(dx^2+dy^2))/(c0*dt+sqrt(dx^2+dy^2)))*(Hz(2,2)-tempz(1,1));
%     % ���Ͻ�
%     Hz(1,Ny+1)=tempz(2,Ny)+((c0*dt-sqrt(dx^2+dy^2))/(c0*dt+sqrt(dx^2+dy^2)))*(Hz(2,Ny)-tempz(1,Ny+1));
%     % ���Ͻ�
%     Hz(Nx+1,Ny+1)=tempz(Nx,Ny)+((c0*dt-sqrt(dx^2+dy^2))/(c0*dt+sqrt(dx^2+dy^2)))*(Hz(Nx,Ny)-tempz(Nx+1,Ny+1));
%     % ���½�
%     Hz(Nx+1,1)=tempz(Nx,2)+((c0*dt-sqrt(dx^2+dy^2))/(c0*dt+sqrt(dx^2+dy^2)))*(Hz(Nx,2)-tempz(Nx+1,1));
    % �ǵ㣨��������
    % �÷�����Ҫ��������ʱ���ǰ�ĳ�ǿ������ר��������һ��ȫ�ֱ���
    % ���½�
    Hz(1,1) = (1-sqrt(2)*c0*dt/sqrt(dx^2+dy^2))*temp2z(1,1)+sqrt(2)*c0*dt/sqrt(dx^2+dy^2)*temp2z(2,2);
    % ���Ͻ�
    Hz(1,Ny+1) = (1-sqrt(2)*c0*dt/sqrt(dx^2+dy^2))*temp2z(1,Ny+1)+sqrt(2)*c0*dt/sqrt(dx^2+dy^2)*temp2z(2,Ny);
    % ���Ͻ�
    Hz(Nx+1,Ny+1) = (1-sqrt(2)*c0*dt/sqrt(dx^2+dy^2))*temp2z(Nx+1,Ny+1)+sqrt(2)*c0*dt/sqrt(dx^2+dy^2)*temp2z(Nx,Ny);
    % ���½�
    Hz(Nx+1,1) = (1-sqrt(2)*c0*dt/sqrt(dx^2+dy^2))*temp2z(Nx+1,1)+sqrt(2)*c0*dt/sqrt(dx^2+dy^2)*temp2z(Nx,2);
    
    temp2z = tempz;
    
    % Visualize fields ���ӻ���    
    imagesc(Hz');
    shading flat;caxis([-1.0 1.0]);axis image;axis xy; 
    title(['Hz, step ',int2str(n)]);xlabel('i');ylabel('j');
    drawnow;

end

%***********************************************************************
% END TIME-STEPPING LOOP ����ѭ��
%***********************************************************************