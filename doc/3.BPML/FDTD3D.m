% This program demonstrates a three-dimensional FDTD simulation.
% The program mainly simulates the propagation of electromagnetic wave in Z 
% directionin free space. The excitation source is a harmonic field source.
% The BPML boundary condition is used as the boundary condition.

% 该程序演示了三维FDTD模拟。
% 该程序主要模拟了电磁波在自由空间中z方向上的传播，激励源为时谐场源，边界条件采用了
% BPML边界条件。

clc;
clear;
close all;

%***********************************************************************
% Fundamental constants 基本常数
%***********************************************************************

eps0 = 8.85e-12;	% permittivity of vacuum 真空介电常数
mu0 = 4*pi*1e-7;	% permeability of vacuum 真空磁导率
c0 = 1/sqrt(mu0*eps0);	% speed of light 光速
z0 = sqrt(mu0/eps0);    % Wave impedance of vacuum 真空中的波阻抗

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

bpml = 8;
m = 4;
sigExmax = (m+1)/(sqrt(epsR)*150*pi*dx);
sigEymax = (m+1)/(sqrt(epsR)*150*pi*dy);
sigEzmax = (m+1)/(sqrt(epsR)*150*pi*dz);
sigHxmax = (m+1)/(sqrt(epsR)*150*pi*dx);
sigHymax = (m+1)/(sqrt(epsR)*150*pi*dy);
sigHzmax = (m+1)/(sqrt(epsR)*150*pi*dz);
% 对于Ez和Hz所需要的分量
sigEzx1 = zeros(Nx+2*bpml,1);
sigEzy1 = zeros(Ny+2*bpml,1);
sigEzx2 = zeros(Nx+2*bpml+1,1);
sigEzy2 = zeros(Ny+2*bpml+1,1);
sigHzx1 = zeros(Nx+2*bpml,1);
sigHzy1 = zeros(Ny+2*bpml,1);
sigHzx2 = zeros(Nx+2*bpml+1,1);
sigHzy2 = zeros(Ny+2*bpml+1,1);
for i = 1:bpml
    sigEzx1(bpml+1-i) = sigExmax*(i/bpml)^m; % 左
    sigEzy1(bpml+1-i) = sigEymax*(i/bpml)^m; % 下
    sigEzx1(Nx+bpml+i) = sigExmax*(i/bpml)^m;	% 右
    sigEzy1(Ny+bpml+i) = sigEymax*(i/bpml)^m;	% 上
    sigEzx2(bpml+1-i) = sigExmax*(i/bpml)^m; % 左
    sigEzy2(bpml+1-i) = sigEymax*(i/bpml)^m; % 下
    sigEzx2(Nx+bpml+1+i) = sigExmax*(i/bpml)^m;	% 右
    sigEzy2(Ny+bpml+1+i) = sigEymax*(i/bpml)^m;	% 上
    sigHzx1(bpml+1-i) = sigHxmax*(i/bpml)^m; % 左
    sigHzy1(bpml+1-i) = sigHymax*(i/bpml)^m; % 下
    sigHzx1(Nx+bpml+i) = sigHxmax*(i/bpml)^m;	% 右
    sigHzy1(Ny+bpml+i) = sigHymax*(i/bpml)^m;	% 上
    sigHzx2(bpml+1-i) = sigHxmax*(i/bpml)^m; % 左
    sigHzy2(bpml+1-i) = sigHymax*(i/bpml)^m; % 下
    sigHzx2(Nx+bpml+1+i) = sigHxmax*(i/bpml)^m;	% 右
    sigHzy2(Ny+bpml+1+i) = sigHymax*(i/bpml)^m;	% 上
end

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

Ex = zeros(Nx+bpml*2,Ny+bpml*2+1,Nz+bpml*2+1);
Ey = zeros(Nx+bpml*2+1,Ny+bpml*2,Nz+bpml*2+1);
Ez = zeros(Nx+bpml*2+1,Ny+bpml*2+1,Nz+bpml*2);
Hx = zeros(Nx+bpml*2+1,Ny+bpml*2,Nz+bpml*2);
Hy = zeros(Nx+bpml*2,Ny+bpml*2+1,Nz+bpml*2);
Hz = zeros(Nx+bpml*2,Ny+bpml*2,Nz+bpml*2+1);
Exy = zeros(Nx+bpml*2,Ny+bpml*2+1,Nz+bpml*2+1);
Exz = zeros(Nx+bpml*2,Ny+bpml*2+1,Nz+bpml*2+1);
Eyx = zeros(Nx+bpml*2+1,Ny+bpml*2,Nz+bpml*2+1);
Eyz = zeros(Nx+bpml*2+1,Ny+bpml*2,Nz+bpml*2+1);
Ezx = zeros(Nx+bpml*2+1,Ny+bpml*2+1,Nz+bpml*2);
Ezy = zeros(Nx+bpml*2+1,Ny+bpml*2+1,Nz+bpml*2);
Hxy = zeros(Nx+bpml*2+1,Ny+bpml*2,Nz+bpml*2);
Hxz = zeros(Nx+bpml*2+1,Ny+bpml*2,Nz+bpml*2);
Hyx = zeros(Nx+bpml*2,Ny+bpml*2+1,Nz+bpml*2);
Hyz = zeros(Nx+bpml*2,Ny+bpml*2+1,Nz+bpml*2);
Hzx = zeros(Nx+bpml*2,Ny+bpml*2,Nz+bpml*2+1);
Hzy = zeros(Nx+bpml*2,Ny+bpml*2,Nz+bpml*2+1);

%***********************************************************************
% BEGIN TIME-STEPPING LOOP 开始循环
%***********************************************************************

for n=1:Nt

    % Set excitation source 设置激励源
    Ez(Jx,Jy,1:Nz) = sin(2*pi*fre*n*dt);
    
    % Update magnetic field 更新磁场
    % Hx
    for i = 1+bpml:Nx+bpml
        for j = 1+bpml:Ny+bpml+1
            for k  = 1+bpml:Nz+bpml+1
                Hx(i,j,k) = CP*Hx(i,j,k)-CQ* ...
                            ((Ez(i,j+1,k)-Ez(i,j,k))/dy- ...
                            (Ey(i,j,k+1)-Ey(i,j,k))/dz);
            end
        end
    end
    % Hy
    for i = 1+bpml:Nx+bpml+1
        for j = 1+bpml:Ny+bpml
            for k  = 1+bpml:Nz+bpml+1
                Hy(i,j,k) = CP*Hy(i,j,k)-CQ* ...
                            ((Ex(i,j,k+1)-Ex(i,j,k))/dz- ...
                            (Ez(i+1,j,k)-Ez(i,j,k))/dx);
            end
        end
    end
    % Hz
    for i = 1+bpml:Nx+bpml+1
        for j = 1+bpml:Ny+bpml+1
            for k  = 1+bpml:Nz+bpml
                Hz(i,j,k) = CP*Hz(i,j,k)-CQ* ...
                            ((Ey(i+1,j,k)-Ey(i,j,k))/dx- ...
                            (Ex(i,j+1,k)-Ex(i,j,k))/dy);
            end
        end
    end
    % Hz分裂成Hzx和Hzy
    for k  = 1+bpml:Nz+bpml+1
%         Hzx(i,j) = exp(-sigHzx2(i)*dt/mu0)*Hzx(i,j)-(1-exp(-sigHzx2(i)*dt/mu0))/(dx*sigHzx2(i))*(Eyz(i,j)-Eyz(i-1,j)+Eyx(i,j)-Eyx(i-1,j));
%         Hzy(i,j) = Hzy(i,j)+1/(2*z0)*(Exy(i,j)-Exy(i,j-1)+Exz(i,j)-Exz(i,j-1));
%         Hzx(i,j) = Hzx(i,j)-1/(2*z0)*(Eyz(i,j)-Eyz(i-1,j)+Eyx(i,j)-Eyx(i-1,j));
%         Hzy(i,j) = exp(-sigHzy2(j)*dt/mu0)*Hzy(i,j)+(1-exp(-sigHzy2(j)*dt/mu0))/(dy*sigHzy2(j))*(Exy(i,j)-Exy(i,j-1)+Exz(i,j)-Exz(i,j-1));
        % 左
        for i = 2:bpml
            for j = 1+bpml:Ny+bpml+1
                Hzx(i,j) = exp(-sigHzx2(i)*dt/mu0)*Hzx(i,j)-(1-exp(-sigHzx2(i)*dt/mu0))/(dx*sigHzx2(i))*(Eyz(i,j)-Eyz(i-1,j)+Eyx(i,j)-Eyx(i-1,j));
                Hzy(i,j) = Hzy(i,j)+1/(2*z0)*(Exy(i,j)-Exy(i,j-1)+Exz(i,j)-Exz(i,j-1));
                Hz(i,j) = Hzx(i,j)+Hzy(i,j);
            end
        end
        % 右
        for i = Nx+bpml+2:Nx+2*bpml
            for j = 1+bpml:Ny+bpml+1
                Hzx(i,j) = exp(-sigHzx2(i)*dt/mu0)*Hzx(i,j)-(1-exp(-sigHzx2(i)*dt/mu0))/(dx*sigHzx2(i))*(Eyz(i,j)-Eyz(i-1,j)+Eyx(i,j)-Eyx(i-1,j));
                Hzy(i,j) = Hzy(i,j)+1/(2*z0)*(Exy(i,j)-Exy(i,j-1)+Exz(i,j)-Exz(i,j-1));
                Hz(i,j) = Hzx(i,j)+Hzy(i,j);
            end
        end
        % 下
        for i = 1+bpml:Nx+bpml+1
            for j = 2:bpml
                Hzx(i,j) = Hzx(i,j)-1/(2*z0)*(Eyz(i,j)-Eyz(i-1,j)+Eyx(i,j)-Eyx(i-1,j));
                Hzy(i,j) = exp(-sigHzy2(j)*dt/mu0)*Hzy(i,j)+(1-exp(-sigHzy2(j)*dt/mu0))/(dy*sigHzy2(j))*(Exy(i,j)-Exy(i,j-1)+Exz(i,j)-Exz(i,j-1));
                Hz(i,j) = Hzx(i,j)+Hzy(i,j);
            end
        end
        % 上
        for i = 1+bpml:Nx+bpml+1
            for j = Ny+bpml+2:Ny+2*bpml
                Hzx(i,j) = Hzx(i,j)-1/(2*z0)*(Eyz(i,j)-Eyz(i-1,j)+Eyx(i,j)-Eyx(i-1,j));
                Hzy(i,j) = exp(-sigHzy2(j)*dt/mu0)*Hzy(i,j)+(1-exp(-sigHzy2(j)*dt/mu0))/(dy*sigHzy2(j))*(Exy(i,j)-Exy(i,j-1)+Exz(i,j)-Exz(i,j-1));
                Hz(i,j) = Hzx(i,j)+Hzy(i,j);
            end
        end
        % 四个角
        % 左下
        for i = 2:bpml
            for j = 2:bpml
                Hzx(i,j) = exp(-sigHzx2(i)*dt/mu0)*Hzx(i,j)-(1-exp(-sigHzx2(i)*dt/mu0))/(dx*sigHzx2(i))*(Eyz(i,j)-Eyz(i-1,j)+Eyx(i,j)-Eyx(i-1,j));
                Hzy(i,j) = exp(-sigHzy2(j)*dt/mu0)*Hzy(i,j)+(1-exp(-sigHzy2(j)*dt/mu0))/(dy*sigHzy2(j))*(Exy(i,j)-Exy(i,j-1)+Exz(i,j)-Exz(i,j-1));
                Hz(i,j) = Hzx(i,j)+Hzy(i,j);
            end
        end
        % 左上
        for i = 2:bpml
            for j = Ny+bpml+2:Ny+2*bpml
                Hzx(i,j) = exp(-sigHzx2(i)*dt/mu0)*Hzx(i,j)-(1-exp(-sigHzx2(i)*dt/mu0))/(dx*sigHzx2(i))*(Eyz(i,j)-Eyz(i-1,j)+Eyx(i,j)-Eyx(i-1,j));
                Hzy(i,j) = exp(-sigHzy2(j)*dt/mu0)*Hzy(i,j)+(1-exp(-sigHzy2(j)*dt/mu0))/(dy*sigHzy2(j))*(Exy(i,j)-Exy(i,j-1)+Exz(i,j)-Exz(i,j-1));
                Hz(i,j) = Hzx(i,j)+Hzy(i,j);
            end
        end
        % 右上
        for i = Nx+bpml+2:Nx+2*bpml
            for j = Ny+bpml+2:Ny+2*bpml
                Hzx(i,j) = exp(-sigHzx2(i)*dt/mu0)*Hzx(i,j)-(1-exp(-sigHzx2(i)*dt/mu0))/(dx*sigHzx2(i))*(Eyz(i,j)-Eyz(i-1,j)+Eyx(i,j)-Eyx(i-1,j));
                Hzy(i,j) = exp(-sigHzy2(j)*dt/mu0)*Hzy(i,j)+(1-exp(-sigHzy2(j)*dt/mu0))/(dy*sigHzy2(j))*(Exy(i,j)-Exy(i,j-1)+Exz(i,j)-Exz(i,j-1));
                Hz(i,j) = Hzx(i,j)+Hzy(i,j);
            end
        end
        % 右下
        for i = Nx+bpml+2:Nx+2*bpml
            for j = 2:bpml
                Hzx(i,j) = exp(-sigHzx2(i)*dt/mu0)*Hzx(i,j)-(1-exp(-sigHzx2(i)*dt/mu0))/(dx*sigHzx2(i))*(Eyz(i,j)-Eyz(i-1,j)+Eyx(i,j)-Eyx(i-1,j));
                Hzy(i,j) = exp(-sigHzy2(j)*dt/mu0)*Hzy(i,j)+(1-exp(-sigHzy2(j)*dt/mu0))/(dy*sigHzy2(j))*(Exy(i,j)-Exy(i,j-1)+Exz(i,j)-Exz(i,j-1));
                Hz(i,j) = Hzx(i,j)+Hzy(i,j);
            end
        end
    end
    
    
    % Update electric field 更新电场
    % Ex
    for i = 1+bpml:Nx+bpml+1
        for j = 1+bpml:Ny+bpml
            for k  = 1+bpml:Nz+bpml
                Ex(i,j,k) = CA*Ex(i,j,k)+CB* ...
                            ((Hz(i,j,k)-Hz(i,j-1,k))/dy- ...
                            (Hy(i,j,k)-Hy(i,j,k-1))/dz);
            end
        end
    end
    % Ey
    for i = 1+bpml:Nx+bpml
        for j = 1+bpml:Ny+bpml+1
            for k  = 1+bpml:Nz+bpml
                Ey(i,j,k) = CA*Ey(i,j,k)+CB* ...
                            ((Hx(i,j,k)-Hx(i,j,k-1))/dz- ...
                            (Hz(i,j,k)-Hz(i-1,j,k))/dx);
            end
        end
    end
    % Ez
    for i = 1+bpml:Nx+bpml
        for j = 1+bpml:Ny+bpml
            for k  = 1+bpml:Nz+bpml+1
                Ez(i,j,k) = CA*Ez(i,j,k)+CB* ...
                            ((Hy(i,j,k)-Hy(i-1,j,k))/dx- ...
                            (Hx(i,j,k)-Hx(i,j-1,k))/dy);
            end
        end
    end
    % Ez分裂成Ezx和Ezy
    for k   = 1+bpml:Nz+bpml+1
%         Ezx(i,j) = exp(-sigEzx1(i)*dt/eps0)*Ezx(i,j)+(1-exp(-sigEzx1(i)*dt/eps0))/(dx*sigEzx1(i))*(Hyz(i,j)-Hyz(i-1,j)+Hyx(i,j)-Hyx(i-1,j));
%         Ezy(i,j) = Ezy(i,j)-1/(2*z0)*(Hxy(i,j)-Hxy(i,j-1)+Hxz(i,j)-Hxz(i,j-1));
%         Ezx(i,j) = Ezx(i,j)+1/(2*z0)*(Hyz(i,j)-Hyz(i-1,j)+Hyx(i,j)-Hyx(i-1,j));
%         Ezy(i,j) = exp(-sigEzy1(j)*dt/eps0)*Ezy(i,j)-(1-exp(-sigEzy1(j)*dt/eps0))/(dy*sigEzy1(j))*(Hxy(i,j)-Hxy(i,j-1)+Hxz(i,j)-Hxz(i,j-1));
%         Ez(i,j) = Ezx(i,j)+Ezy(i,j);
        % 左
        for i = 2:bpml
            for j = 1+bpml:Ny+bpml+1
                Ezx(i,j) = exp(-sigEzx1(i)*dt/eps0)*Ezx(i,j)+(1-exp(-sigEzx1(i)*dt/eps0))/(dx*sigEzx1(i))*(Hyz(i,j)-Hyz(i-1,j)+Hyx(i,j)-Hyx(i-1,j));
                Ezy(i,j) = Ezy(i,j)-1/(2*z0)*(Hxy(i,j)-Hxy(i,j-1)+Hxz(i,j)-Hxz(i,j-1));
                Ez(i,j) = Ezx(i,j)+Ezy(i,j);
            end
        end
        % 右
        for i = Nx+bpml+2:Nx+2*bpml
            for j = 1+bpml:Ny+bpml+1
                Ezx(i,j) = exp(-sigEzx1(i)*dt/eps0)*Ezx(i,j)+(1-exp(-sigEzx1(i)*dt/eps0))/(dx*sigEzx1(i))*(Hyz(i,j)-Hyz(i-1,j)+Hyx(i,j)-Hyx(i-1,j));
                Ezy(i,j) = Ezy(i,j)-1/(2*z0)*(Hxy(i,j)-Hxy(i,j-1)+Hxz(i,j)-Hxz(i,j-1));
                Ez(i,j) = Ezx(i,j)+Ezy(i,j);
            end
        end
        % 下
        for i = 1+bpml:Nx+bpml+1
            for j = 2:bpml
                Ezx(i,j) = Ezx(i,j)+1/(2*z0)*(Hyz(i,j)-Hyz(i-1,j)+Hyx(i,j)-Hyx(i-1,j));
                Ezy(i,j) = exp(-sigEzy1(j)*dt/eps0)*Ezy(i,j)-(1-exp(-sigEzy1(j)*dt/eps0))/(dy*sigEzy1(j))*(Hxy(i,j)-Hxy(i,j-1)+Hxz(i,j)-Hxz(i,j-1));
                Ez(i,j) = Ezx(i,j)+Ezy(i,j);
            end
        end
        % 上
        for i = 1+bpml:Nx+bpml+1
            for j = Ny+bpml+2:Ny+2*bpml
                Ezx(i,j) = Ezx(i,j)+1/(2*z0)*(Hyz(i,j)-Hyz(i-1,j)+Hyx(i,j)-Hyx(i-1,j));
                Ezy(i,j) = exp(-sigEzy1(j)*dt/eps0)*Ezy(i,j)-(1-exp(-sigEzy1(j)*dt/eps0))/(dy*sigEzy1(j))*(Hxy(i,j)-Hxy(i,j-1)+Hxz(i,j)-Hxz(i,j-1));
                Ez(i,j) = Ezx(i,j)+Ezy(i,j);
            end
        end
        % 四个角
        % 左下
        for i = 2:bpml
            for j = 2:bpml
                Ezx(i,j) = exp(-sigEzx1(i)*dt/eps0)*Ezx(i,j)+(1-exp(-sigEzx1(i)*dt/eps0))/(dx*sigEzx1(i))*(Hyz(i,j)-Hyz(i-1,j)+Hyx(i,j)-Hyx(i-1,j));
                Ezy(i,j) = exp(-sigEzy1(j)*dt/eps0)*Ezy(i,j)-(1-exp(-sigEzy1(j)*dt/eps0))/(dy*sigEzy1(j))*(Hxy(i,j)-Hxy(i,j-1)+Hxz(i,j)-Hxz(i,j-1));
                Ez(i,j) = Ezx(i,j)+Ezy(i,j);
            end
        end
        % 左上
        for i = 2:bpml
            for j = Ny+bpml+2:Ny+2*bpml
                Ezx(i,j) = exp(-sigEzx1(i)*dt/eps0)*Ezx(i,j)+(1-exp(-sigEzx1(i)*dt/eps0))/(dx*sigEzx1(i))*(Hyz(i,j)-Hyz(i-1,j)+Hyx(i,j)-Hyx(i-1,j));
                Ezy(i,j) = exp(-sigEzy1(j)*dt/eps0)*Ezy(i,j)-(1-exp(-sigEzy1(j)*dt/eps0))/(dy*sigEzy1(j))*(Hxy(i,j)-Hxy(i,j-1)+Hxz(i,j)-Hxz(i,j-1));
                Ez(i,j) = Ezx(i,j)+Ezy(i,j);
            end
        end
        % 右上
        for i = Nx+bpml+2:Nx+2*bpml
            for j = Ny+bpml+2:Ny+2*bpml
                Ezx(i,j) = exp(-sigEzx1(i)*dt/eps0)*Ezx(i,j)+(1-exp(-sigEzx1(i)*dt/eps0))/(dx*sigEzx1(i))*(Hyz(i,j)-Hyz(i-1,j)+Hyx(i,j)-Hyx(i-1,j));
                Ezy(i,j) = exp(-sigEzy1(j)*dt/eps0)*Ezy(i,j)-(1-exp(-sigEzy1(j)*dt/eps0))/(dy*sigEzy1(j))*(Hxy(i,j)-Hxy(i,j-1)+Hxz(i,j)-Hxz(i,j-1));
                Ez(i,j) = Ezx(i,j)+Ezy(i,j);
            end
        end
        % 右下
        for i = Nx+bpml+2:Nx+2*bpml
            for j = 2:bpml
                Ezx(i,j) = exp(-sigEzx1(i)*dt/eps0)*Ezx(i,j)+(1-exp(-sigEzx1(i)*dt/eps0))/(dx*sigEzx1(i))*(Hyz(i,j)-Hyz(i-1,j)+Hyx(i,j)-Hyx(i-1,j));
                Ezy(i,j) = exp(-sigEzy1(j)*dt/eps0)*Ezy(i,j)-(1-exp(-sigEzy1(j)*dt/eps0))/(dy*sigEzy1(j))*(Hxy(i,j)-Hxy(i,j-1)+Hxz(i,j)-Hxz(i,j-1));
                Ez(i,j) = Ezx(i,j)+Ezy(i,j);
            end
        end
    end
         
    % Set boundary conditions 设置边界条件

    % Visualize fields 可视化场
    imagesc(Ez(:,:,Jz)');
    shading flat;caxis([-1.0 1.0]);axis image;axis xy; 
    title(['Ez(i,j,k=',int2str(Jz),'),step = ',int2str(n)]);xlabel('i');ylabel('j');
    pause(0.01);

end

%***********************************************************************
% END TIME-STEPPING LOOP 结束循环
%***********************************************************************