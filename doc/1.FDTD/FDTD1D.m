% This program demonstrates a one-dimensional FDTD simulation.
% The program mainly simulates the propagation of electromagnetic wave in Z 
% directionin free space. The excitation source is a harmonic field source.
% The simplest PEC boundary condition is used as the boundary condition.

% 该程序演示了一维FDTD模拟。该程序主要模拟电磁波在自由空间中z方向上的传播。
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

Nz = 100;	% number of cells in 1D problem space 一维问题空间中的单元数
Nt = 150;	% number of iterations 迭代次数
dz = 3e-2;	% space step 空间步长
dt = 1/(2.0*c0*sqrt(1/dz^2));	%time step 时间步长

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
Jz = round(Nz/2);   % position 位置

%***********************************************************************
% Initializing field arrays 初始化场
%***********************************************************************

Ex = zeros(Nz+1,1);
Hy = zeros(Nz,1);

%***********************************************************************
% BEGIN TIME-STEPPING LOOP 开始迭代
%***********************************************************************

for n=1:Nt
    
    % Set excitation source 设置激励源
%     Ex(1) = sin(2.0*pi*fre*dt*n);
%     Ex(Jz) = 1;
    Ex(Jz) = sin(2.0*pi*fre*dt*n);
%     Ex(Nz+1) = sin(2.0*pi*fre*dt*n);

    % Update magnetic field 更新磁场
    for k = 1:Nz
        Hy(k) = CP*Hy(k)-CQ*(Ex(k+1)-Ex(k))/dz;
    end
    
    % Update electric field 更新电场
    for k = 2:Nz
        Ex(k) = CA*Ex(k)-CB*(Hy(k)-Hy(k-1))/dz;
    end
   
    % Set boundary conditions 设置边界条件
    
    % Visualize fields 可视化场
    subplot(2,1,1),plot(Ex,'r');
    xlabel('z');ylabel(['Ex, step ',int2str(n)]);axis([0 100 -2 2]);
    subplot(2,1,2),plot(Hy,'b');
    xlabel('z');ylabel(['Hy, step ',int2str(n)]);axis([0 100 -5e-3 5e-3]);
    drawnow;

end

%***********************************************************************
% END TIME-STEPPING LOOP 终止迭代
%***********************************************************************