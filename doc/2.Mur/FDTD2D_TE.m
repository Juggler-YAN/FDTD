% This program demonstrates a two-dimensional FDTD simulation(Modified TE).
% The program mainly simulates the propagation of electromagnetic wave in 
% free space. The excitation source is a harmonic field source. The Mur
% absorbing boundary condition is used as the boundary condition.

% 该程序演示了二维FDTD模拟（修改后的TE）。该程序主要模拟电磁波在自由空间中
% 的传播。激励源为时谐场源。边界条件为Mur吸收边界条件。

clc;
clear;
close all;

%***********************************************************************
% Fundamental constants 基本常数
%***********************************************************************

eps0 = 8.85e-12;	% permittivity of vacuum 真空介电常数
mu0 = 4*pi*1e-7;	% permeability of vacuum 真空磁导率
c0 = 1/sqrt(mu0*eps0);	% speed of light 光速

%***********************************************************************
% Mesh parameters 网格参数
%***********************************************************************

Nx = 50;	% number of cells in 2D problem space 二维问题空间中的单元数
Ny = 50;
Nt = 150;	% number of iterations 迭代次数
dx = 3e-2;	% space step 空间步长
dy = 3e-2;
dt = 1/(2.0*c0*sqrt(1/dx^2+1/dy^2));	%time step 时间步长

%***********************************************************************
% Material properties 媒质特性
%***********************************************************************

epsR = 1;	% relative permittivity 相对介电常数
muR = 1;	% relative permeability 相对电导率
sigE = 0;	% electric conductivity 电导率
sigH = 0;	% effective magnetism conductivity 等效磁导率

%***********************************************************************
% Boundary conditions 边界条件
%***********************************************************************

% PEC边界条件

%***********************************************************************
% Updating coefficients 更新系数
%***********************************************************************

CA = (eps0*epsR/dt-sigE/2)/(eps0*epsR/dt+sigE/2);
CB = 1/(eps0*epsR/dt+sigE/2);
CP = (mu0*muR/dt-sigH/2)/(mu0*muR/dt+sigH/2);
CQ = 1/(mu0*muR/dt+sigH/2);

%***********************************************************************
% Source excitation 源激励
%***********************************************************************

fre = 1.0e+9; % frequency 频率
Jx = round(Nx/2);	% position 位置
Jy = round(Ny/2);

%***********************************************************************
% Initializing field arrays 初始化场
%***********************************************************************

Ex = zeros(Nx+1,Ny);
Ey = zeros(Nx,Ny+1);
Hz = zeros(Nx+1,Ny+1);

%***********************************************************************
% BEGIN TIME-STEPPING LOOP 开始迭代
%***********************************************************************

temp2z = Hz;
for n=1:Nt
    
    % Set excitation source 设置激励源
    Hz(Jx,Jy) = sin(2*pi*fre*n*dt);
    
    % Update electric field 更新磁场
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
    
    % Update magnetic field 更新电场
    % Hz
    for i = 2:Nx
        for j = 2:Ny
           Hz(i,j) = CP*Hz(i,j)-CQ* ...
                      ((Ey(i,j)-Ey(i-1,j))/dx-(Ex(i,j)-Ex(i,j-1))/dy);
        end
    end 
    
    % Set boundary conditions 设置边界条件
    % 边
    % 一阶Mur
%     % 左
%     for j = 2:Ny
%         Hz(1,j) = tempz(2,j)+((c0*dt-dx)/(c0*dt+dx))*(Hz(2,j)-tempz(1,j));
%     end
%     % 右
%     for j = 2:Ny
%         Hz(Nx+1,j) = tempz(Nx,j)+((c0*dt-dx)/(c0*dt+dx))*(Hz(Nx,j)-tempz(Nx+1,j));
%     end
%     % 下
%     for i = 2:Nx
%         Hz(i,1) = tempz(i,2)+((c0*dt-dy)/(c0*dt+dy))*(Hz(i,2)-tempz(i,1));
%     end
%     % 上
%     for i = 2:Nx
%         Hz(i,Ny+1) = tempz(i,Ny)+((c0*dt-dy)/(c0*dt+dy))*(Hz(i,Ny)-tempz(i,Ny+1));
%     end
    % 二阶Mur
    % 左
    for j = 2:Ny
        Hz(1,j) = tempz(2,j)+((c0*dt-dx)/(c0*dt+dx))*(Hz(2,j)-tempz(1,j))+ ...
                  (c0*c0*eps0*epsR*dt)/(2*(c0*dt+dx))*(dx/dy)*(Ex(1,j)-Ex(1,j-1)+Ex(2,j)-Ex(2,j-1));
    end
    % 右
    for j = 2:Ny
        Hz(Nx+1,j) = tempz(Nx,j)+((c0*dt-dx)/(c0*dt+dx))*(Hz(Nx,j)-tempz(Nx+1,j))+ ...
                     (c0*c0*eps0*epsR*dt)/(2*(c0*dt+dx))*(dx/dy)*(Ex(Nx+1,j)-Ex(Nx+1,j-1)+Ex(Nx,j)-Ex(Nx,j-1));
    end
    % 下
    for i = 2:Nx
        Hz(i,1) = tempz(i,2)+((c0*dt-dy)/(c0*dt+dy))*(Hz(i,2)-tempz(i,1))- ...
                  (c0*c0*eps0*epsR*dt)/(2*(c0*dt+dy))*(dy/dx)*(Ey(i,1)-Ey(i-1,1)+Ey(i,2)-Ey(i-1,2));
    end
    % 上
    for i = 2:Nx
        Hz(i,Ny+1) = tempz(i,Ny)+((c0*dt-dy)/(c0*dt+dy))*(Hz(i,Ny)-tempz(i,Ny+1))- ...
                     (c0*c0*eps0*epsR*dt)/(2*(c0*dt+dy))*(dy/dx)*(Ey(i,Ny+1)-Ey(i-1,Ny+1)+Ey(i,Ny)-Ey(i-1,Ny));
    end
    % 角点
	% 方法一
%     % 左下
%     Hz(1,1)=tempz(2,2)+((c0*dt-sqrt(dx^2+dy^2))/(c0*dt+sqrt(dx^2+dy^2)))*(Hz(2,2)-tempz(1,1));
%     % 左上
%     Hz(1,Ny+1)=tempz(2,Ny)+((c0*dt-sqrt(dx^2+dy^2))/(c0*dt+sqrt(dx^2+dy^2)))*(Hz(2,Ny)-tempz(1,Ny+1));
%     % 右上
%     Hz(Nx+1,Ny+1)=tempz(Nx,Ny)+((c0*dt-sqrt(dx^2+dy^2))/(c0*dt+sqrt(dx^2+dy^2)))*(Hz(Nx,Ny)-tempz(Nx+1,Ny+1));
%     % 右下
%     Hz(Nx+1,1)=tempz(Nx,2)+((c0*dt-sqrt(dx^2+dy^2))/(c0*dt+sqrt(dx^2+dy^2)))*(Hz(Nx,2)-tempz(Nx+1,1));
    % 方法二
    % 左下
    Hz(1,1) = (1-sqrt(2)*c0*dt/sqrt(dx^2+dy^2))*temp2z(1,1)+sqrt(2)*c0*dt/sqrt(dx^2+dy^2)*temp2z(2,2);
    % 左上
    Hz(1,Ny+1) = (1-sqrt(2)*c0*dt/sqrt(dx^2+dy^2))*temp2z(1,Ny+1)+sqrt(2)*c0*dt/sqrt(dx^2+dy^2)*temp2z(2,Ny);
    % 右上
    Hz(Nx+1,Ny+1) = (1-sqrt(2)*c0*dt/sqrt(dx^2+dy^2))*temp2z(Nx+1,Ny+1)+sqrt(2)*c0*dt/sqrt(dx^2+dy^2)*temp2z(Nx,Ny);
    % 右下
    Hz(Nx+1,1) = (1-sqrt(2)*c0*dt/sqrt(dx^2+dy^2))*temp2z(Nx+1,1)+sqrt(2)*c0*dt/sqrt(dx^2+dy^2)*temp2z(Nx,2);
    
    temp2z = tempz;
    
    % Visualize fields 可视化场 
    imagesc(Hz');
    shading flat;caxis([-1.0 1.0]);axis image;axis xy; 
    title(['Hz, step ',int2str(n)]);xlabel('i');ylabel('j');
    drawnow;

end

%***********************************************************************
% END TIME-STEPPING LOOP 终止循环
%***********************************************************************