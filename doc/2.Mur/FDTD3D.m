% This program demonstrates a three-dimensional FDTD simulation.
% The program mainly simulates the propagation of electromagnetic wave in 
% free space. The excitation source is a harmonic field source.
% The Mur absorbing boundary condition is used as the boundary condition.

% 该程序演示了三维FDTD模拟。该程序主要模拟电磁波在自由空间中的传播�?
% �?��源为时谐场源。边界条件为Mur吸收边界条件�?

clc;
clear;
close all;

%***********************************************************************
% Fundamental constants 基本常数
%***********************************************************************

eps0 = 8.85e-12;	% permittivity of vacuum 真空介电常数
mu0 = 4*pi*1e-7;	% permeability of vacuum 真空磁导�?
c0 = 1/sqrt(mu0*eps0);	% speed of light 光�?

%***********************************************************************
% Mesh parameters 网格参数
%***********************************************************************

Nx = 50;	% number of cells in 3D problem space 三维问题空间中的单元�?
Ny = 50;
Nz = 50;
Nt = 150;	% number of iterations 迭代次数
dx = 3e-2;	% space step 空间步长
dy = 3e-2;
dz = 3e-2;
dt = 1/(2.0*c0*sqrt(1/dx^2+1/dy^2+1/dz^2));	%time step 时间步长

%***********************************************************************
% Material properties 媒质特�?
%***********************************************************************

epsR = 1;	% relative permittivity 相对介电常数
muR = 1;	% relative permeability 相对磁导�?
sigE = 0;	% electric conductivity 电导�?
sigH = 0;	% effective magnetism conductivity 等效磁导�?

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
% Source excitation �?���?
%***********************************************************************

fre = 1.0e+9; % frequency 频率
Jx = round(Nx/2);   % position 位置
Jy = round(Ny/2);
Jz = round(Nz/2);

%***********************************************************************
% Initializing field arrays 初始化场
%***********************************************************************

Ex = zeros(Nx,Ny+1,Nz+1);
Ey = zeros(Nx+1,Ny,Nz+1);
Ez = zeros(Nx+1,Ny+1,Nz);
Hx = zeros(Nx+1,Ny,Nz);
Hy = zeros(Nx,Ny+1,Nz);
Hz = zeros(Nx,Ny,Nz+1);

%***********************************************************************
% BEGIN TIME-STEPPING LOOP �?��迭代
%***********************************************************************

temp2x = Ex;
temp2y = Ey;
temp2z = Ez;
for n=1:Nt

    % Set excitation source 设置�?���?
    Ez(Jx,Jy,1:Nz) = sin(2*pi*fre*n*dt);
    
    % Update magnetic field 更新磁场
    % Hx
    for i = 2:Nx
        for j = 1:Ny
            for k  = 1:Nz
                Hx(i,j,k) = CP*Hx(i,j,k)-CQ* ...
                            ((Ez(i,j+1,k)-Ez(i,j,k))/dy- ...
                            (Ey(i,j,k+1)-Ey(i,j,k))/dz);
            end
        end
    end
    % Hy
    for i = 1:Nx
        for j = 2:Ny
            for k  = 1:Nz
                Hy(i,j,k) = CP*Hy(i,j,k)-CQ* ...
                            ((Ex(i,j,k+1)-Ex(i,j,k))/dz- ...
                            (Ez(i+1,j,k)-Ez(i,j,k))/dx);
            end
        end
    end
    % Hz
    for i = 1:Nx
        for j = 1:Ny
            for k  = 2:Nz
                Hz(i,j,k) = CP*Hz(i,j,k)-CQ* ...
                            ((Ey(i+1,j,k)-Ey(i,j,k))/dx- ...
                            (Ex(i,j+1,k)-Ex(i,j,k))/dy);
            end
        end
    end
    
    tempx = Ex;
    tempy = Ey;
    tempz = Ez;
    
    % Update electric field 更新电场
    % Ex
    for i = 1:Nx
        for j = 2:Ny
            for k  = 2:Nz
                Ex(i,j,k) = CA*Ex(i,j,k)+CB* ...
                            ((Hz(i,j,k)-Hz(i,j-1,k))/dy- ...
                            (Hy(i,j,k)-Hy(i,j,k-1))/dz);
            end
        end
    end
    % Ey
    for i = 2:Nx
        for j = 1:Ny
            for k  = 2:Nz
                Ey(i,j,k) = CA*Ey(i,j,k)+CB* ...
                            ((Hx(i,j,k)-Hx(i,j,k-1))/dz- ...
                            (Hz(i,j,k)-Hz(i-1,j,k))/dx);
            end
        end
    end
    % Ez
    for i = 2:Nx
        for j = 2:Ny
            for k  = 1:Nz
                Ez(i,j,k) = CA*Ez(i,j,k)+CB* ...
                            ((Hy(i,j,k)-Hy(i-1,j,k))/dx- ...
                            (Hx(i,j,k)-Hx(i,j-1,k))/dy);
            end
        end
    end
         
    % Set boundary conditions 设置边界条件
    % �?��Mur
    % �?
    % �?
    for j = 2:Ny-1
        for k = 2:Nz
            Ey(Nx+1,j,k) = tempy(Nx,j,k)+((c0*dt-dx)/(c0*dt+dx))*(Ey(Nx,j,k)-tempy(Nx+1,j,k));
        end
    end
    for j = 2:Ny
        for k = 2:Nz-1
            Ez(Nx+1,j,k) = tempz(Nx,j,k)+((c0*dt-dx)/(c0*dt+dx))*(Ez(Nx,j,k)-tempz(Nx+1,j,k));
        end
    end
    % �?
    for j = 2:Ny-1
        for k = 2:Nz
            Ey(1,j,k) = tempy(2,j,k)+((c0*dt-dx)/(c0*dt+dx))*(Ey(2,j,k)-tempy(1,j,k));
        end
    end
    for j = 2:Ny
        for k = 2:Nz-1
            Ez(1,j,k) = tempz(2,j,k)+((c0*dt-dx)/(c0*dt+dx))*(Ez(2,j,k)-tempz(1,j,k));
        end
    end
    % �?
    for i = 2:Nx-1
        for k = 2:Nz
            Ex(i,1,k) = tempx(i,2,k)+((c0*dt-dy)/(c0*dt+dy))*(Ex(i,2,k)-tempx(i,1,k));
        end
    end
    for i = 2:Nx
        for k = 2:Nz-1
            Ez(i,1,k) = tempz(i,2,k)+((c0*dt-dy)/(c0*dt+dy))*(Ez(i,2,k)-tempz(i,1,k));
        end
    end
    % �?
    for i = 2:Nx-1
        for k = 2:Nz
            Ex(i,Ny+1,k) = tempx(i,Ny,k)+((c0*dt-dy)/(c0*dt+dy))*(Ex(i,Ny,k)-tempx(i,Ny+1,k));
        end
    end
    for i = 2:Nx
        for k = 2:Nz-1
            Ez(i,Ny+1,k) = tempz(i,Ny,k)+((c0*dt-dy)/(c0*dt+dy))*(Ez(i,Ny,k)-tempz(i,Ny+1,k));
        end
    end
    % �?
    for i = 2:Nx-1
        for j = 2:Ny
            Ex(i,j,Nz+1) = tempx(i,j,Nz)+((c0*dt-dz)/(c0*dt+dz))*(Ex(i,j,Nz)-tempx(i,j,Nz+1));
        end
    end
    for i = 2:Nx
        for j = 2:Ny-1
            Ey(i,j,Nz+1) = tempy(i,j,Nz)+((c0*dt-dz)/(c0*dt+dz))*(Ey(i,j,Nz)-tempy(i,j,Nz+1));
        end
    end
    % �?
    for i = 2:Nx-1
        for j = 2:Ny
            Ex(i,j,1) = tempx(i,j,2)+((c0*dt-dz)/(c0*dt+dz))*(Ex(i,j,2)-tempx(i,j,1));
        end
    end
    for i = 2:Nx
        for j = 2:Ny-1
            Ey(i,j,1) = tempy(i,j,2)+((c0*dt-dz)/(c0*dt+dz))*(Ey(i,j,2)-tempy(i,j,1));
        end
    end
    % �?
    % x方向
    for i = 2:Nx-1
        Ex(i,1,1) = tempx(i,2,2)+((c0*dt-sqrt(dy^2+dz^2))/(c0*dt+sqrt(dy^2+dz^2)))*(Ex(i,2,2)-tempx(i,1,1));
        Ex(i,1,Nz+1) = tempx(i,2,Nz)+((c0*dt-sqrt(dy^2+dz^2))/(c0*dt+sqrt(dy^2+dz^2)))*(Ex(i,2,Nz)-tempx(i,1,Nz+1));
        Ex(i,Ny+1,Nz+1) = tempx(i,Ny,Nz)+((c0*dt-sqrt(dy^2+dz^2))/(c0*dt+sqrt(dy^2+dz^2)))*(Ex(i,Ny,Nz)-tempx(i,Ny+1,Nz+1));
        Ex(i,Ny+1,1) = tempx(i,Ny,2)+((c0*dt-sqrt(dy^2+dz^2))/(c0*dt+sqrt(dy^2+dz^2)))*(Ex(i,Ny,2)-tempx(i,Ny+1,1));
    end
    % y方向
    for j = 2:Ny-1
        Ey(1,j,1) = tempy(2,j,2)+((c0*dt-sqrt(dx^2+dz^2))/(c0*dt+sqrt(dx^2+dz^2)))*(Ey(2,j,2)-tempy(1,j,1));
        Ey(1,j,Nz+1) = tempy(2,j,Nz)+((c0*dt-sqrt(dx^2+dz^2))/(c0*dt+sqrt(dx^2+dz^2)))*(Ey(2,j,Nz)-tempy(1,j,Nz+1));
        Ey(Nx+1,j,Nz+1) = tempy(Nx,j,Nz)+((c0*dt-sqrt(dx^2+dz^2))/(c0*dt+sqrt(dx^2+dz^2)))*(Ey(Nx,j,Nz)-tempy(Nx+1,j,Nz+1));
        Ey(Nx+1,j,1) = tempy(Nx,j,2)+((c0*dt-sqrt(dx^2+dz^2))/(c0*dt+sqrt(dx^2+dz^2)))*(Ey(Nx,j,2)-tempy(Nx+1,j,1));
    end
    % z方向
    for k = 2:Nz-1
        Ez(1,1,k) = tempz(2,2,k)+((c0*dt-sqrt(dx^2+dy^2))/(c0*dt+sqrt(dx^2+dy^2)))*(Ez(2,2,k)-tempz(1,1,k));
        Ez(1,Ny+1,k) = tempz(2,Ny,k)+((c0*dt-sqrt(dx^2+dy^2))/(c0*dt+sqrt(dx^2+dy^2)))*(Ez(2,Ny,k)-tempz(1,Ny+1,k));
        Ez(Nx+1,Ny+1,k) = tempz(Nx,Ny,k)+((c0*dt-sqrt(dx^2+dy^2))/(c0*dt+sqrt(dx^2+dy^2)))*(Ez(Nx,Ny,k)-tempz(Nx+1,Ny+1,k));
        Ez(Nx+1,1,k) = tempz(Nx,2,k)+((c0*dt-sqrt(dx^2+dy^2))/(c0*dt+sqrt(dx^2+dy^2)))*(Ez(Nx,2,k)-tempz(Nx+1,1,k));
    end
    % �?
    % 角顶点处无电磁场节点分量，无�?���?

    % 二阶Mur
%     % �?
%     % �?
%     for j = 2:Ny-1
%         for k = 2:Nz
%             Ey(Nx+1,j,k) = -temp2y(Nx,j,k)+((c0*dt-dx)/(c0*dt+dx))*(Ey(Nx,j,k)+temp2y(Nx+1,j,k))+ ...
%                            2*dx/(c0*dt+dx)*(tempy(Nx+1,j,k)+tempy(Nx,j,k))+ ...
%                            (dx*(c0*dt)^2)/(2*(dy^2)*(c0*dt+dx))* ...
%                            (tempy(Nx+1,j+1,k)-2*tempy(Nx+1,j,k)+tempy(Nx+1,j-1,k)+tempy(Nx,j+1,k)-2*tempy(Nx,j,k)+tempy(Nx,j-1,k))+ ...
%                            (dx*(c0*dt)^2)/(2*(dz^2)*(c0*dt+dx))* ...
%                            (tempy(Nx+1,j,k+1)-2*tempy(Nx+1,j,k)+tempy(Nx+1,j,k-1)+tempy(Nx,j,k+1)-2*tempy(Nx,j,k)+tempy(Nx,j,k-1));
%         end
%     end
%     for j = 2:Ny
%         for k = 2:Nz-1
%             Ez(Nx+1,j,k) = -temp2z(Nx,j,k)+((c0*dt-dx)/(c0*dt+dx))*(Ez(Nx,j,k)+temp2z(Nx+1,j,k))+ ...
%                            2*dx/(c0*dt+dx)*(tempz(Nx+1,j,k)+tempz(Nx,j,k))+ ...
%                            (dx*(c0*dt)^2)/(2*(dy^2)*(c0*dt+dx))* ...
%                            (tempz(Nx+1,j+1,k)-2*tempz(Nx+1,j,k)+tempz(Nx+1,j-1,k)+tempz(Nx,j+1,k)-2*tempz(Nx,j,k)+tempz(Nx,j-1,k))+ ...
%                            (dx*(c0*dt)^2)/(2*(dz^2)*(c0*dt+dx))* ...
%                            (tempz(Nx+1,j,k+1)-2*tempz(Nx+1,j,k)+tempz(Nx+1,j,k-1)+tempz(Nx,j,k+1)-2*tempz(Nx,j,k)+tempz(Nx,j,k-1));
%                         
%         end
%     end
%     % �?
%     for j = 2:Ny-1
%         for k = 2:Nz
%             Ey(1,j,k) = -temp2y(2,j,k)+((c0*dt-dx)/(c0*dt+dx))*(Ey(2,j,k)+temp2y(1,j,k))+ ...
%                         2*dx/(c0*dt+dx)*(tempy(1,j,k)+tempy(2,j,k))+ ...
%                         (dx*(c0*dt)^2)/(2*(dy^2)*(c0*dt+dx))* ...
%                         (tempy(1,j+1,k)-2*tempy(1,j,k)+tempy(1,j-1,k)+tempy(2,j+1,k)-2*tempy(2,j,k)+tempy(2,j-1,k))+ ...
%                         (dx*(c0*dt)^2)/(2*(dz^2)*(c0*dt+dx))* ...
%                         (tempy(1,j,k+1)-2*tempy(1,j,k)+tempy(1,j,k-1)+tempy(2,j,k+1)-2*tempy(2,j,k)+tempy(2,j,k-1));
%         end
%     end
%     for j = 2:Ny
%         for k = 2:Nz-1
%             Ez(1,j,k) = -temp2z(2,j,k)+((c0*dt-dx)/(c0*dt+dx))*(Ez(2,j,k)+temp2z(1,j,k))+ ...
%                         2*dx/(c0*dt+dx)*(tempz(1,j,k)+tempz(2,j,k))+ ...
%                         (dx*(c0*dt)^2)/(2*(dy^2)*(c0*dt+dx))* ...
%                         (tempz(1,j+1,k)-2*tempz(1,j,k)+tempz(1,j-1,k)+tempz(2,j+1,k)-2*tempz(2,j,k)+tempz(2,j-1,k))+ ...
%                         (dx*(c0*dt)^2)/(2*(dz^2)*(c0*dt+dx))* ...
%                         (tempz(1,j,k+1)-2*tempz(1,j,k)+tempz(1,j,k-1)+tempz(2,j,k+1)-2*tempz(2,j,k)+tempz(2,j,k-1));
%         end
%     end
%     % �?
%     for i = 2:Nx-1
%         for k = 2:Nz
%             Ex(i,1,k) = -temp2x(i,2,k)+((c0*dt-dy)/(c0*dt+dy))*(Ex(i,2,k)+temp2x(i,1,k))+ ...
%                         2*dy/(c0*dt+dy)*(tempx(i,1,k)+tempx(i,2,k))+ ...
%                         (dy*(c0*dt)^2)/(2*(dx^2)*(c0*dt+dy))* ...
%                         (tempx(i+1,1,k)-2*tempx(i,1,k)+tempx(i-1,1,k)+tempx(i+1,2,k)-2*tempx(i,2,k)+tempx(i-1,2,k))+ ...
%                         (dy*(c0*dt)^2)/(2*(dz^2)*(c0*dt+dy))* ...
%                         (tempx(i,1,k+1)-2*tempx(i,1,k)+tempx(i,1,k-1)+tempx(i,2,k+1)-2*tempx(i,2,k)+tempx(i,2,k-1));
%         end
%     end
%     for i = 2:Nx
%         for k = 2:Nz-1
%             Ez(i,1,k) = -temp2z(i,2,k)+((c0*dt-dy)/(c0*dt+dy))*(Ez(i,2,k)+temp2z(i,1,k))+ ...
%                         2*dy/(c0*dt+dy)*(tempz(i,1,k)+tempz(i,2,k))+ ...
%                         (dy*(c0*dt)^2)/(2*(dx^2)*(c0*dt+dy))* ...
%                         (tempz(i+1,1,k)-2*tempz(i,1,k)+tempz(i-1,1,k)+tempz(i+1,2,k)-2*tempz(i,2,k)+tempz(i-1,2,k))+ ...
%                         (dy*(c0*dt)^2)/(2*(dz^2)*(c0*dt+dy))* ...
%                         (tempz(i,1,k+1)-2*tempz(i,1,k)+tempz(i,1,k-1)+tempz(i,2,k+1)-2*tempz(i,2,k)+tempz(i,2,k-1));
%         end
%     end
%     % �?
%     for i = 2:Nx-1
%         for k = 2:Nz
%             Ex(i,Ny+1,k) = -temp2x(i,Ny,k)+((c0*dt-dy)/(c0*dt+dy))*(Ex(i,Ny,k)+temp2x(i,Ny+1,k))+ ...
%                            2*dy/(c0*dt+dy)*(tempx(i,Ny+1,k)+tempx(i,Ny,k))+ ...
%                            (dy*(c0*dt)^2)/(2*(dx^2)*(c0*dt+dy))* ...
%                            (tempx(i+1,Ny+1,k)-2*tempx(i,Ny+1,k)+tempx(i-1,Ny+1,k)+tempx(i+1,Ny,k)-2*tempx(i,Ny,k)+tempx(i-1,Ny,k))+ ...
%                            (dy*(c0*dt)^2)/(2*(dz^2)*(c0*dt+dy))* ...
%                            (tempx(i,Ny+1,k+1)-2*tempx(i,Ny+1,k)+tempx(i,Ny+1,k-1)+tempx(i,Ny,k+1)-2*tempx(i,Ny,k)+tempx(i,Ny,k-1));
%         end
%     end
%     for i = 2:Nx
%         for k = 2:Nz-1
%             Ez(i,Ny+1,k) = -temp2z(i,Ny,k)+((c0*dt-dy)/(c0*dt+dy))*(Ez(i,Ny,k)+temp2z(i,Ny+1,k))+ ...
%                            2*dy/(c0*dt+dy)*(tempz(i,Ny+1,k)+tempz(i,Ny,k))+ ...
%                            (dy*(c0*dt)^2)/(2*(dx^2)*(c0*dt+dy))* ...
%                            (tempz(i+1,Ny+1,k)-2*tempz(i,Ny+1,k)+tempz(i-1,Ny+1,k)+tempz(i+1,Ny,k)-2*tempz(i,Ny,k)+tempz(i-1,Ny,k))+ ...
%                            (dy*(c0*dt)^2)/(2*(dz^2)*(c0*dt+dy))* ...
%                            (tempz(i,Ny+1,k+1)-2*tempz(i,Ny+1,k)+tempz(i,Ny+1,k-1)+tempz(i,Ny,k+1)-2*tempz(i,Ny,k)+tempz(i,Ny,k-1));
%         end
%     end
%     % �?
%     for i = 2:Nx-1
%         for j = 2:Ny
%             Ex(i,j,Nz+1) = -temp2x(i,j,Nz)+((c0*dt-dz)/(c0*dt+dz))*(Ex(i,j,Nz)+temp2x(i,j,Nz+1))+ ...
%                            2*dz/(c0*dt+dz)*(tempx(i,j,Nz+1)+tempx(i,j,Nz))+ ...
%                            (dz*(c0*dt)^2)/(2*(dx^2)*(c0*dt+dz))* ...
%                            (tempx(i+1,j,Nz+1)-2*tempx(i,j,Nz+1)+tempx(i-1,j,Nz+1)+tempx(i+1,j,Nz)-2*tempx(i,j,Nz)+tempx(i-1,j,Nz))+ ...
%                            (dz*(c0*dt)^2)/(2*(dy^2)*(c0*dt+dz))* ...
%                            (tempx(i,j+1,Nz+1)-2*tempx(i,j,Nz+1)+tempx(i,j-1,Nz+1)+tempx(i,j+1,Nz)-2*tempx(i,j,Nz)+tempx(i,j-1,Nz));
%         end
%     end
%     for i = 2:Nx
%         for j = 2:Ny-1
%             Ey(i,j,Nz+1) = -temp2y(i,j,Nz)+((c0*dt-dz)/(c0*dt+dz))*(Ey(i,j,Nz)+temp2y(i,j,Nz+1))+ ...
%                            2*dz/(c0*dt+dz)*(tempy(i,j,Nz+1)+tempy(i,j,Nz))+ ...
%                            (dz*(c0*dt)^2)/(2*(dx^2)*(c0*dt+dz))* ...
%                            (tempy(i+1,j,Nz+1)-2*tempy(i,j,Nz+1)+tempy(i-1,j,Nz+1)+tempy(i+1,j,Nz)-2*tempy(i,j,Nz)+tempy(i-1,j,Nz))+ ...
%                            (dz*(c0*dt)^2)/(2*(dy^2)*(c0*dt+dz))* ...
%                            (tempy(i,j+1,Nz+1)-2*tempy(i,j,Nz+1)+tempy(i,j-1,Nz+1)+tempy(i,j+1,Nz)-2*tempy(i,j,Nz)+tempy(i,j-1,Nz));
%         end
%     end
%     % �?
%     for i = 2:Nx-1
%         for j = 2:Ny
%             Ex(i,j,1) = -temp2x(i,j,2)+((c0*dt-dz)/(c0*dt+dz))*(Ex(i,j,2)+temp2x(i,j,1))+ ...
%                         2*dz/(c0*dt+dz)*(tempx(i,j,1)+tempx(i,j,2))+ ...
%                         (dz*(c0*dt)^2)/(2*(dx^2)*(c0*dt+dz))* ...
%                         (tempx(i+1,j,1)-2*tempx(i,j,1)+tempx(i-1,j,1)+tempx(i+1,j,2)-2*tempx(i,j,2)+tempx(i-1,j,2))+ ...
%                         (dz*(c0*dt)^2)/(2*(dy^2)*(c0*dt+dz))* ...
%                         (tempx(i,j+1,1)-2*tempx(i,j,1)+tempx(i,j-1,1)+tempx(i,j+1,2)-2*tempx(i,j,2)+tempx(i,j-1,2));
%         end
%     end
%     for i = 2:Nx
%         for j = 2:Ny-1
%             Ey(i,j,1) = -temp2y(i,j,2)+((c0*dt-dz)/(c0*dt+dz))*(Ey(i,j,2)+temp2y(i,j,1))+ ...
%                         2*dz/(c0*dt+dz)*(tempy(i,j,1)+tempy(i,j,2))+ ...
%                         (dz*(c0*dt)^2)/(2*(dx^2)*(c0*dt+dz))* ...
%                         (tempy(i+1,j,1)-2*tempy(i,j,1)+tempy(i-1,j,1)+tempy(i+1,j,2)-2*tempy(i,j,2)+tempy(i-1,j,2))+ ...
%                         (dz*(c0*dt)^2)/(2*(dy^2)*(c0*dt+dz))* ...
%                         (tempy(i,j+1,1)-2*tempy(i,j,1)+tempy(i,j-1,1)+tempy(i,j+1,2)-2*tempy(i,j,2)+tempy(i,j-1,2));
%         end
%     end
%     % �?
%     % x方向
%     for i = 2:Nx-1
%         Ex(i,1,1) = tempx(i,2,2)+((c0*dt-sqrt(dy^2+dz^2))/(c0*dt+sqrt(dy^2+dz^2)))*(Ex(i,2,2)-tempx(i,1,1));
%         Ex(i,1,Nz+1) = tempx(i,2,Nz)+((c0*dt-sqrt(dy^2+dz^2))/(c0*dt+sqrt(dy^2+dz^2)))*(Ex(i,2,Nz)-tempx(i,1,Nz+1));
%         Ex(i,Ny+1,Nz+1) = tempx(i,Ny,Nz)+((c0*dt-sqrt(dy^2+dz^2))/(c0*dt+sqrt(dy^2+dz^2)))*(Ex(i,Ny,Nz)-tempx(i,Ny+1,Nz+1));
%         Ex(i,Ny+1,1) = tempx(i,Ny,2)+((c0*dt-sqrt(dy^2+dz^2))/(c0*dt+sqrt(dy^2+dz^2)))*(Ex(i,Ny,2)-tempx(i,Ny+1,1));
%     end
%     % y方向
%     for j = 2:Ny-1
%         Ey(1,j,1) = tempy(2,j,2)+((c0*dt-sqrt(dx^2+dz^2))/(c0*dt+sqrt(dx^2+dz^2)))*(Ey(2,j,2)-tempy(1,j,1));
%         Ey(1,j,Nz+1) = tempy(2,j,Nz)+((c0*dt-sqrt(dx^2+dz^2))/(c0*dt+norm([dx,dz])))*(Ey(2,j,Nz)-tempy(1,j,Nz+1));
%         Ey(Nx+1,j,Nz+1) = tempy(Nx,j,Nz)+((c0*dt-sqrt(dx^2+dz^2))/(c0*dt+sqrt(dx^2+dz^2)))*(Ey(Nx,j,Nz)-tempy(Nx+1,j,Nz+1));
%         Ey(Nx+1,j,1) = tempy(Nx,j,2)+((c0*dt-sqrt(dx^2+dz^2))/(c0*dt+sqrt(dx^2+dz^2)))*(Ey(Nx,j,2)-tempy(Nx+1,j,1));
%     end
%     % z方向
%     for k = 2:Nz-1
%         Ez(1,1,k) = tempz(2,2,k)+((c0*dt-sqrt(dx^2+dy^2))/(c0*dt+sqrt(dx^2+dy^2)))*(Ez(2,2,k)-tempz(1,1,k));
%         Ez(1,Ny+1,k) = tempz(2,Ny,k)+((c0*dt-sqrt(dx^2+dy^2))/(c0*dt+sqrt(dx^2+dy^2)))*(Ez(2,Ny,k)-tempz(1,Ny+1,k));
%         Ez(Nx+1,Ny+1,k) = tempz(Nx,Ny,k)+((c0*dt-sqrt(dx^2+dy^2))/(c0*dt+sqrt(dx^2+dy^2)))*(Ez(Nx,Ny,k)-tempz(Nx+1,Ny+1,k));
%         Ez(Nx+1,1,k) = tempz(Nx,2,k)+((c0*dt-sqrt(dx^2+dy^2))/(c0*dt+sqrt(dx^2+dy^2)))*(Ez(Nx,2,k)-tempz(Nx+1,1,k));
%     end
%     % �?
%     % 角顶点处无电磁场节点分量，无�?���?

    % 二阶Mur优化
%     % �?
%     % �?
%     for j = 3:Ny-2
%         for k = 3:Nz-1
%             Ey(Nx+1,j,k) = -temp2y(Nx,j,k)+((c0*dt-dx)/(c0*dt+dx))*(Ey(Nx,j,k)+temp2y(Nx+1,j,k))+ ...
%                            2*dx/(c0*dt+dx)*(tempy(Nx+1,j,k)+tempy(Nx,j,k))+ ...
%                            (dx*(c0*dt)^2)/(2*(dy^2)*(c0*dt+dx))* ...
%                            (tempy(Nx+1,j+1,k)-2*tempy(Nx+1,j,k)+tempy(Nx+1,j-1,k)+tempy(Nx,j+1,k)-2*tempy(Nx,j,k)+tempy(Nx,j-1,k))+ ...
%                            (dx*(c0*dt)^2)/(2*(dz^2)*(c0*dt+dx))* ...
%                            (tempy(Nx+1,j,k+1)-2*tempy(Nx+1,j,k)+tempy(Nx+1,j,k-1)+tempy(Nx,j,k+1)-2*tempy(Nx,j,k)+tempy(Nx,j,k-1));
%         end
%     end
%     for j = 3:Ny-1
%         for k = 3:Nz-2
%             Ez(Nx+1,j,k) = -temp2z(Nx,j,k)+((c0*dt-dx)/(c0*dt+dx))*(Ez(Nx,j,k)+temp2z(Nx+1,j,k))+ ...
%                            2*dx/(c0*dt+dx)*(tempz(Nx+1,j,k)+tempz(Nx,j,k))+ ...
%                            (dx*(c0*dt)^2)/(2*(dy^2)*(c0*dt+dx))* ...
%                            (tempz(Nx+1,j+1,k)-2*tempz(Nx+1,j,k)+tempz(Nx+1,j-1,k)+tempz(Nx,j+1,k)-2*tempz(Nx,j,k)+tempz(Nx,j-1,k))+ ...
%                            (dx*(c0*dt)^2)/(2*(dz^2)*(c0*dt+dx))* ...
%                            (tempz(Nx+1,j,k+1)-2*tempz(Nx+1,j,k)+tempz(Nx+1,j,k-1)+tempz(Nx,j,k+1)-2*tempz(Nx,j,k)+tempz(Nx,j,k-1));
%                         
%         end
%     end
%     % �?
%     for j = 3:Ny-2
%         for k = 3:Nz-1
%             Ey(1,j,k) = -temp2y(2,j,k)+((c0*dt-dx)/(c0*dt+dx))*(Ey(2,j,k)+temp2y(1,j,k))+ ...
%                         2*dx/(c0*dt+dx)*(tempy(1,j,k)+tempy(2,j,k))+ ...
%                         (dx*(c0*dt)^2)/(2*(dy^2)*(c0*dt+dx))* ...
%                         (tempy(1,j+1,k)-2*tempy(1,j,k)+tempy(1,j-1,k)+tempy(2,j+1,k)-2*tempy(2,j,k)+tempy(2,j-1,k))+ ...
%                         (dx*(c0*dt)^2)/(2*(dz^2)*(c0*dt+dx))* ...
%                         (tempy(1,j,k+1)-2*tempy(1,j,k)+tempy(1,j,k-1)+tempy(2,j,k+1)-2*tempy(2,j,k)+tempy(2,j,k-1));
%         end
%     end
%     for j = 3:Ny-1
%         for k = 3:Nz-2
%             Ez(1,j,k) = -temp2z(2,j,k)+((c0*dt-dx)/(c0*dt+dx))*(Ez(2,j,k)+temp2z(1,j,k))+ ...
%                         2*dx/(c0*dt+dx)*(tempz(1,j,k)+tempz(2,j,k))+ ...
%                         (dx*(c0*dt)^2)/(2*(dy^2)*(c0*dt+dx))* ...
%                         (tempz(1,j+1,k)-2*tempz(1,j,k)+tempz(1,j-1,k)+tempz(2,j+1,k)-2*tempz(2,j,k)+tempz(2,j-1,k))+ ...
%                         (dx*(c0*dt)^2)/(2*(dz^2)*(c0*dt+dx))* ...
%                         (tempz(1,j,k+1)-2*tempz(1,j,k)+tempz(1,j,k-1)+tempz(2,j,k+1)-2*tempz(2,j,k)+tempz(2,j,k-1));
%         end
%     end
%     % �?
%     for i = 3:Nx-2
%         for k = 3:Nz-1
%             Ex(i,1,k) = -temp2x(i,2,k)+((c0*dt-dy)/(c0*dt+dy))*(Ex(i,2,k)+temp2x(i,1,k))+ ...
%                         2*dy/(c0*dt+dy)*(tempx(i,1,k)+tempx(i,2,k))+ ...
%                         (dy*(c0*dt)^2)/(2*(dx^2)*(c0*dt+dy))* ...
%                         (tempx(i+1,1,k)-2*tempx(i,1,k)+tempx(i-1,1,k)+tempx(i+1,2,k)-2*tempx(i,2,k)+tempx(i-1,2,k))+ ...
%                         (dy*(c0*dt)^2)/(2*(dz^2)*(c0*dt+dy))* ...
%                         (tempx(i,1,k+1)-2*tempx(i,1,k)+tempx(i,1,k-1)+tempx(i,2,k+1)-2*tempx(i,2,k)+tempx(i,2,k-1));
%         end
%     end
%     for i = 3:Nx-1
%         for k = 3:Nz-2
%             Ez(i,1,k) = -temp2z(i,2,k)+((c0*dt-dy)/(c0*dt+dy))*(Ez(i,2,k)+temp2z(i,1,k))+ ...
%                         2*dy/(c0*dt+dy)*(tempz(i,1,k)+tempz(i,2,k))+ ...
%                         (dy*(c0*dt)^2)/(2*(dx^2)*(c0*dt+dy))* ...
%                         (tempz(i+1,1,k)-2*tempz(i,1,k)+tempz(i-1,1,k)+tempz(i+1,2,k)-2*tempz(i,2,k)+tempz(i-1,2,k))+ ...
%                         (dy*(c0*dt)^2)/(2*(dz^2)*(c0*dt+dy))* ...
%                         (tempz(i,1,k+1)-2*tempz(i,1,k)+tempz(i,1,k-1)+tempz(i,2,k+1)-2*tempz(i,2,k)+tempz(i,2,k-1));
%         end
%     end
%     % �?
%     for i = 3:Nx-2
%         for k = 3:Nz-1
%             Ex(i,Ny+1,k) = -temp2x(i,Ny,k)+((c0*dt-dy)/(c0*dt+dy))*(Ex(i,Ny,k)+temp2x(i,Ny+1,k))+ ...
%                            2*dy/(c0*dt+dy)*(tempx(i,Ny+1,k)+tempx(i,Ny,k))+ ...
%                            (dy*(c0*dt)^2)/(2*(dx^2)*(c0*dt+dy))* ...
%                            (tempx(i+1,Ny+1,k)-2*tempx(i,Ny+1,k)+tempx(i-1,Ny+1,k)+tempx(i+1,Ny,k)-2*tempx(i,Ny,k)+tempx(i-1,Ny,k))+ ...
%                            (dy*(c0*dt)^2)/(2*(dz^2)*(c0*dt+dy))* ...
%                            (tempx(i,Ny+1,k+1)-2*tempx(i,Ny+1,k)+tempx(i,Ny+1,k-1)+tempx(i,Ny,k+1)-2*tempx(i,Ny,k)+tempx(i,Ny,k-1));
%         end
%     end
%     for i = 3:Nx-1
%         for k = 3:Nz-2
%             Ez(i,Ny+1,k) = -temp2z(i,Ny,k)+((c0*dt-dy)/(c0*dt+dy))*(Ez(i,Ny,k)+temp2z(i,Ny+1,k))+ ...
%                            2*dy/(c0*dt+dy)*(tempz(i,Ny+1,k)+tempz(i,Ny,k))+ ...
%                            (dy*(c0*dt)^2)/(2*(dx^2)*(c0*dt+dy))* ...
%                            (tempz(i+1,Ny+1,k)-2*tempz(i,Ny+1,k)+tempz(i-1,Ny+1,k)+tempz(i+1,Ny,k)-2*tempz(i,Ny,k)+tempz(i-1,Ny,k))+ ...
%                            (dy*(c0*dt)^2)/(2*(dz^2)*(c0*dt+dy))* ...
%                            (tempz(i,Ny+1,k+1)-2*tempz(i,Ny+1,k)+tempz(i,Ny+1,k-1)+tempz(i,Ny,k+1)-2*tempz(i,Ny,k)+tempz(i,Ny,k-1));
%         end
%     end
%     % �?
%     for i = 3:Nx-2
%         for j = 3:Ny-1
%             Ex(i,j,Nz+1) = -temp2x(i,j,Nz)+((c0*dt-dz)/(c0*dt+dz))*(Ex(i,j,Nz)+temp2x(i,j,Nz+1))+ ...
%                            2*dz/(c0*dt+dz)*(tempx(i,j,Nz+1)+tempx(i,j,Nz))+ ...
%                            (dz*(c0*dt)^2)/(2*(dx^2)*(c0*dt+dz))* ...
%                            (tempx(i+1,j,Nz+1)-2*tempx(i,j,Nz+1)+tempx(i-1,j,Nz+1)+tempx(i+1,j,Nz)-2*tempx(i,j,Nz)+tempx(i-1,j,Nz))+ ...
%                            (dz*(c0*dt)^2)/(2*(dy^2)*(c0*dt+dz))* ...
%                            (tempx(i,j+1,Nz+1)-2*tempx(i,j,Nz+1)+tempx(i,j-1,Nz+1)+tempx(i,j+1,Nz)-2*tempx(i,j,Nz)+tempx(i,j-1,Nz));
%         end
%     end
%     for i = 3:Nx-1
%         for j = 3:Ny-2
%             Ey(i,j,Nz+1) = -temp2y(i,j,Nz)+((c0*dt-dz)/(c0*dt+dz))*(Ey(i,j,Nz)+temp2y(i,j,Nz+1))+ ...
%                            2*dz/(c0*dt+dz)*(tempy(i,j,Nz+1)+tempy(i,j,Nz))+ ...
%                            (dz*(c0*dt)^2)/(2*(dx^2)*(c0*dt+dz))* ...
%                            (tempy(i+1,j,Nz+1)-2*tempy(i,j,Nz+1)+tempy(i-1,j,Nz+1)+tempy(i+1,j,Nz)-2*tempy(i,j,Nz)+tempy(i-1,j,Nz))+ ...
%                            (dz*(c0*dt)^2)/(2*(dy^2)*(c0*dt+dz))* ...
%                            (tempy(i,j+1,Nz+1)-2*tempy(i,j,Nz+1)+tempy(i,j-1,Nz+1)+tempy(i,j+1,Nz)-2*tempy(i,j,Nz)+tempy(i,j-1,Nz));
%         end
%     end
%     % �?
%     for i = 3:Nx-2
%         for j = 3:Ny-1
%             Ex(i,j,1) = -temp2x(i,j,2)+((c0*dt-dz)/(c0*dt+dz))*(Ex(i,j,2)+temp2x(i,j,1))+ ...
%                         2*dz/(c0*dt+dz)*(tempx(i,j,1)+tempx(i,j,2))+ ...
%                         (dz*(c0*dt)^2)/(2*(dx^2)*(c0*dt+dz))* ...
%                         (tempx(i+1,j,1)-2*tempx(i,j,1)+tempx(i-1,j,1)+tempx(i+1,j,2)-2*tempx(i,j,2)+tempx(i-1,j,2))+ ...
%                         (dz*(c0*dt)^2)/(2*(dy^2)*(c0*dt+dz))* ...
%                         (tempx(i,j+1,1)-2*tempx(i,j,1)+tempx(i,j-1,1)+tempx(i,j+1,2)-2*tempx(i,j,2)+tempx(i,j-1,2));
%         end
%     end
%     for i = 3:Nx-1
%         for j = 3:Ny-2
%             Ey(i,j,1) = -temp2y(i,j,2)+((c0*dt-dz)/(c0*dt+dz))*(Ey(i,j,2)+temp2y(i,j,1))+ ...
%                         2*dz/(c0*dt+dz)*(tempy(i,j,1)+tempy(i,j,2))+ ...
%                         (dz*(c0*dt)^2)/(2*(dx^2)*(c0*dt+dz))* ...
%                         (tempy(i+1,j,1)-2*tempy(i,j,1)+tempy(i-1,j,1)+tempy(i+1,j,2)-2*tempy(i,j,2)+tempy(i-1,j,2))+ ...
%                         (dz*(c0*dt)^2)/(2*(dy^2)*(c0*dt+dz))* ...
%                         (tempy(i,j+1,1)-2*tempy(i,j,1)+tempy(i,j-1,1)+tempy(i,j+1,2)-2*tempy(i,j,2)+tempy(i,j-1,2));
%         end
%     end
%     % �?
%     % x方向
%     for i = 3:Nx-2
%         Ex(i,2,2) = tempx(i,3,3)+((c0*dt-norm([dy,dz]))/(c0*dt+norm([dy,dz])))*(Ex(i,3,3)-tempx(i,2,2));
%         Ex(i,2,Nz) = tempx(i,3,Nz-1)+((c0*dt-norm([dy,dz]))/(c0*dt+norm([dy,dz])))*(Ex(i,3,Nz-1)-tempx(i,2,Nz));
%         Ex(i,Ny,Nz) = tempx(i,Ny-1,Nz-1)+((c0*dt-norm([dy,dz]))/(c0*dt+norm([dy,dz])))*(Ex(i,Ny-1,Nz-1)-tempx(i,Ny,Nz));
%         Ex(i,Ny,2) = tempx(i,Ny-1,3)+((c0*dt-norm([dy,dz]))/(c0*dt+norm([dy,dz])))*(Ex(i,Ny-1,3)-tempx(i,Ny,2));
%     end
%     % y方向
%     for j = 3:Ny-2
%         Ey(2,j,2) = tempy(3,j,3)+((c0*dt-norm([dx,dz]))/(c0*dt+norm([dx,dz])))*(Ey(3,j,3)-tempy(2,j,2));
%         Ey(2,j,Nz) = tempy(3,j,Nz-1)+((c0*dt-norm([dx,dz]))/(c0*dt+norm([dx,dz])))*(Ey(3,j,Nz-1)-tempy(2,j,Nz));
%         Ey(Nx,j,Nz) = tempy(Nx-1,j,Nz-1)+((c0*dt-norm([dx,dz]))/(c0*dt+norm([dx,dz])))*(Ey(Nx-1,j,Nz-1)-tempy(Nx,j,Nz));
%         Ey(Nx,j,2) = tempy(Nx-1,j,3)+((c0*dt-norm([dx,dz]))/(c0*dt+norm([dx,dz])))*(Ey(Nx-1,j,3)-tempy(Nx,j,2));
%     end
%     % z方向
%     for k = 3:Nz-2
%         Ez(2,2,k) = tempz(3,3,k)+((c0*dt-norm([dx,dy]))/(c0*dt+norm([dx,dy])))*(Ez(3,3,k)-tempz(2,2,k));
%         Ez(2,Ny,k) = tempz(3,Ny-1,k)+((c0*dt-norm([dx,dy]))/(c0*dt+norm([dx,dy])))*(Ez(3,Ny-1,k)-tempz(2,Ny,k));
%         Ez(Nx,Ny,k) = tempz(Nx-1,Ny-1,k)+((c0*dt-norm([dx,dy]))/(c0*dt+norm([dx,dy])))*(Ez(Nx-1,Ny-1,k)-tempz(Nx,Ny,k));
%         Ez(Nx,2,k) = tempz(Nx-1,3,k)+((c0*dt-norm([dx,dy]))/(c0*dt+norm([dx,dy])))*(Ez(Nx-1,3,k)-tempz(Nx,2,k));
%     end
%     % �?
%     % 角顶点处无电磁场节点分量，无�?���?
    
    temp2x = tempx;
    temp2y = tempy;
    temp2z = tempz;

    % Visualize fields 可视化场
    imagesc(Ez(:,:,Jz)');
    shading flat;caxis([-1.0 1.0]);axis image;axis xy; 
    title(['Ez(i,j,k=',int2str(Jz),'),step = ',int2str(n)]);xlabel('i');ylabel('j');
    drawnow;

end

%***********************************************************************
% END TIME-STEPPING LOOP 终止迭代
%***********************************************************************