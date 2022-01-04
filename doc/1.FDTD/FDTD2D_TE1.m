% This program demonstrates a two-dimensional FDTD simulation(TE).
% The program mainly simulates the propagation of electromagnetic wave in 
% free space. The excitation source is a harmonic field source. The
% simplest PEC boundary condition is used as the boundary condition.

% 该程序演示了二维FDTD模拟（TE）。该程序主要模拟电磁波在自由空间中的传播。
% 激励源为时谐场源。边界条件为最简单的PEC边界条件。

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

Ex = zeros(Nx,Ny+1);
Ey = zeros(Nx+1,Ny);
Hz = zeros(Nx,Ny);

%***********************************************************************
% BEGIN TIME-STEPPING LOOP 开始迭代
%***********************************************************************

for n=1:Nt
    
    % Set excitation source 设置激励源
    Hz(Jx,Jy) = sin(2*pi*fre*n*dt);
    
    % Update electric field 更新磁场
    % Ex
    for i = 1:Nx
        for j = 2:Ny
            Ex(i,j) = CA*Ex(i,j)+CB*(Hz(i,j)-Hz(i,j-1))/dy;
        end
    end
    % Ey
    for i = 2:Nx
        for j = 1:Ny
            Ey(i,j) = CA*Ey(i,j)-CB*(Hz(i,j)-Hz(i-1,j))/dx;
        end
    end

    % Update magnetic field 更新电场
    % Hz
    for i = 1:Nx
        for j = 1:Ny
           Hz(i,j) = CP*Hz(i,j)-CQ* ...
                      ((Ey(i+1,j)-Ey(i,j))/dx-(Ex(i,j+1)-Ex(i,j))/dy);
        end
    end
         
    % Set boundary conditions 设置边界条件
    
    % Visualize fields 可视化场
    imagesc(Hz');
    shading flat;caxis([-1.0 1.0]);axis image;axis xy; 
    title(['Hz, step ',int2str(n)]);xlabel('i');ylabel('j');
    drawnow;

end

%***********************************************************************
% END TIME-STEPPING LOOP 终止迭代
%***********************************************************************