% This program demonstrates a three-dimensional FDTD simulation.
% The program mainly simulates the propagation of electromagnetic wave in Z 
% directionin free space. The excitation source is a harmonic field source.
% The simplest PEC boundary condition is used as the boundary condition.

% 该程序演示了三维FDTD模拟。
% 该程序主要模拟了电磁波在自由空间中z方向上的传播，激励源为时谐场源，边界条件
% 采用了最简单的PEC边界条件。

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

Nx = 50;	% number of cells in 3D problem space 三维问题空间中的单元数
Ny = 50;
Nz = 50;
Nt = 150;	% number of iterations 迭代次数
dx = 3e-2;	% space step 空间步长
dy = 3e-2;
dz = 3e-2;
dt = 1/(2.0*c0*sqrt(1/dx^2+1/dy^2+1/dz^2));	%time step 时间步长

%***********************************************************************
% Material properties 媒质特性
%***********************************************************************

epsR = 1;	% relative permittivity 相对介电常数
muR = 1;	% relative permeability 相对磁导率
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

fre = 1.0e+9; % frequency 激励源频率
Jx = round(Nx/2);   % position 激励源位置
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
% BEGIN TIME-STEPPING LOOP 开始循环
%***********************************************************************

for n=1:Nt

    % Set excitation source 设置激励源
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

    % Visualize fields 可视化场
    imagesc(Ez(:,:,Jz)');
    shading flat;caxis([-1.0 1.0]);axis image;axis xy; 
    title(['Ez(i,j,k=',int2str(Jz),'),step = ',int2str(n)]);xlabel('i');ylabel('j');
    pause(0.01)

end

%***********************************************************************
% END TIME-STEPPING LOOP 结束循环
%***********************************************************************