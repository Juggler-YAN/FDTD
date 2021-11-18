% This program demonstrates a two-dimensional FDTD simulation(TM).
% The program mainly simulates the propagation of electromagnetic wave in 
% free space. The excitation source is a harmonic field source. The
% BPML boundary condition is used as the boundary condition.

% 该程序演示了二维FDTD模拟(TM)。
% 该程序主要模拟了电磁波在自由空间中的传播，激励源为时谐场源，边界条件采用了
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

Nx = 50;	% number of cells in 2D problem space 二维问题空间中的单元数
Ny = 50;
Nt = 300;	% number of iterations 迭代次数
dx = 3e-2;	% space step 空间步长
dy = 3e-2;
dt = 1/(2.0*c0*sqrt(1/dx^2+1/dy^2));	%time step 时间步长

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
sigHxmax = (m+1)/(sqrt(epsR)*150*pi*dx);
sigHymax = (m+1)/(sqrt(epsR)*150*pi*dy);
sigEx = zeros(Nx+2*bpml+1,1);
sigEy = zeros(Ny+2*bpml+1,1);
sigHx = zeros(Nx+2*bpml,1);
sigHy = zeros(Ny+2*bpml,1);

for i = 1:bpml
    sigEx(bpml+1-i) = sigExmax*((i*dx)/(bpml*dx))^m; % 左
    sigEy(bpml+1-i) = sigEymax*((i*dy)/(bpml*dy))^m; % 下
    sigEx(Nx+bpml+1+i) = sigExmax*((i*dx)/(bpml*dx))^m;	% 右
    sigEy(Ny+bpml+1+i) = sigEymax*((i*dy)/(bpml*dy))^m;	% 上
    sigHx(bpml+1-i) = sigHxmax*((i*dx)/(bpml*dx))^m; % 左
    sigHy(bpml+1-i) = sigHymax*((i*dy)/(bpml*dy))^m; % 下
    sigHx(Nx+bpml+i) = sigHxmax*((i*dx)/(bpml*dx))^m;	% 右
    sigHy(Ny+bpml+i) = sigHymax*((i*dy)/(bpml*dy))^m;	% 上
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
Jx = round(Nx/2);	% position 激励源位置
Jy = round(Ny/2);

%***********************************************************************
% Initializing field arrays 初始化场
%***********************************************************************

Hx = zeros(Nx+2*bpml+1,Ny+2*bpml);
Hy = zeros(Nx+2*bpml,Ny+2*bpml+1);
Ez = zeros(Nx+2*bpml+1,Ny+2*bpml+1);
Ezx = zeros(Nx+2*bpml+1,Ny+2*bpml+1);
Ezy = zeros(Nx+2*bpml+1,Ny+2*bpml+1);

%***********************************************************************
% BEGIN TIME-STEPPING LOOP 开始循环
%***********************************************************************

for n=1:Nt
    
    % Set excitation source 设置激励源
    Ez(Jx,Jy) = sin(2*pi*fre*n*dt);
    
    % Update magnetic field 更新磁场
    % Hx
    for i = 1+bpml:Nx+bpml+1
        for j = 1+bpml:Ny+bpml
            Hx(i,j) = CP*Hx(i,j)-CQ*(Ez(i,j+1)-Ez(i,j))/dy;
        end
    end
    % Hy
    for i = 1+bpml:Nx+bpml
        for j = 1+bpml:Ny+bpml+1
            Hy(i,j) = CP*Hy(i,j)+CQ*(Ez(i+1,j)-Ez(i,j))/dx;
        end
    end
    % PML
    % 左
    for i = 1:bpml
        for j = 1+bpml:Ny+bpml
            Hx(i,j) = Hx(i,j)-z0/2*(Ez(i,j+1)-Ez(i,j));
        end
    end
    for i = 1:bpml
        for j = 1+bpml:Ny+bpml+1
            Hy(i,j) = exp(-sigHx(i)*dt/mu0)*Hy(i,j)+(1-exp(-sigHx(i)*dt/mu0))/(dx*sigHx(i))*(Ez(i+1,j)-Ez(i,j));
        end
    end
    % 右
    for i = Nx+bpml+2:Nx+2*bpml+1
        for j = 1+bpml:Ny+bpml
            Hx(i,j) = Hx(i,j)-z0/2*(Ez(i,j+1)-Ez(i,j));
        end
    end
    for i = Nx+bpml+1:Nx+2*bpml
        for j = 1+bpml:Ny+bpml+1
            Hy(i,j) = exp(-sigHx(i)*dt/mu0)*Hy(i,j)+(1-exp(-sigHx(i)*dt/mu0))/(dx*sigHx(i))*(Ez(i+1,j)-Ez(i,j));
        end
    end
    % 下
    for i = 1+bpml:Nx+bpml+1
        for j = 1:bpml
            Hx(i,j) = exp(-sigHy(j)*dt/mu0)*Hx(i,j)-(1-exp(-sigHy(j)*dt/mu0))/(dy*sigHy(j))*(Ez(i,j+1)-Ez(i,j));
        end
    end
    for i = 1+bpml:Nx+bpml
        for j = 1:bpml
            Hy(i,j) = Hy(i,j)+z0/2*(Ez(i+1,j)-Ez(i,j));
        end
    end
    % 上
    for i = 1+bpml:Nx+bpml+1
        for j = Ny+bpml+1:Ny+2*bpml
            Hx(i,j) = exp(-sigHy(j)*dt/mu0)*Hx(i,j)-(1-exp(-sigHy(j)*dt/mu0))/(dy*sigHy(j))*(Ez(i,j+1)-Ez(i,j));
        end
    end
    for i = 1+bpml:Nx+bpml
        for j = Ny+bpml+2:Ny+2*bpml
            Hy(i,j) = Hy(i,j)+z0/2*(Ez(i+1,j)-Ez(i,j));
        end
    end
    % 四个角
    % 左下
    for i = 1:bpml
        for j = 1:bpml
            Hx(i,j) = exp(-sigHy(j)*dt/mu0)*Hx(i,j)-(1-exp(-sigHy(j)*dt/mu0))/(dy*sigHy(j))*(Ez(i,j+1)-Ez(i,j));
        end
    end
    for i = 1:bpml
        for j = 1:bpml
            Hy(i,j) = exp(-sigHx(i)*dt/mu0)*Hy(i,j)+(1-exp(-sigHx(i)*dt/mu0))/(dx*sigHx(i))*(Ez(i+1,j)-Ez(i,j));
        end
    end
    % 左上
    for i = 1:bpml
        for j = Ny+bpml+1:Ny+2*bpml
            Hx(i,j) = exp(-sigHy(j)*dt/mu0)*Hx(i,j)-(1-exp(-sigHy(j)*dt/mu0))/(dy*sigHy(j))*(Ez(i,j+1)-Ez(i,j));
        end
    end
    for i = 1:bpml
        for j = Ny+bpml+2:Ny+2*bpml+1
            Hy(i,j) = exp(-sigHx(i)*dt/mu0)*Hy(i,j)+(1-exp(-sigHx(i)*dt/mu0))/(dx*sigHx(i))*(Ez(i+1,j)-Ez(i,j));
        end
    end
    % 右上
    for i = Nx+bpml+2:Nx+2*bpml+1
        for j = Ny+bpml+1:Ny+2*bpml
            Hx(i,j) = exp(-sigHy(j)*dt/mu0)*Hx(i,j)-(1-exp(-sigHy(j)*dt/mu0))/(dy*sigHy(j))*(Ez(i,j+1)-Ez(i,j));
        end
    end
    for i = Nx+bpml+1:Nx+2*bpml
        for j = Ny+bpml+2:Ny+2*bpml+1
            Hy(i,j) = exp(-sigHx(i)*dt/mu0)*Hy(i,j)+(1-exp(-sigHx(i)*dt/mu0))/(dx*sigHx(i))*(Ez(i+1,j)-Ez(i,j));
        end
    end
    % 右下
    for i = Nx+bpml+2:Nx+2*bpml+1
        for j = 1:bpml
            Hx(i,j) = exp(-sigHy(j)*dt/mu0)*Hx(i,j)-(1-exp(-sigHy(j)*dt/mu0))/(dy*sigHy(j))*(Ez(i,j+1)-Ez(i,j));
        end
    end
    for i = Nx+bpml+1:Nx+2*bpml
        for j = 1:bpml
            Hy(i,j) = exp(-sigHx(i)*dt/mu0)*Hy(i,j)+(1-exp(-sigHx(i)*dt/mu0))/(dx*sigHx(i))*(Ez(i+1,j)-Ez(i,j));
        end
    end

%             Hx(i,j) = exp(-sigHy(j)*dt/mu0)*Hx(i,j)-(1-exp(-sigHy(j)*dt/mu0))/(dy*sigHy(j))*(Ez(i,j+1)-Ez(i,j));
%             Hy(i,j) = exp(-sigHx(i)*dt/mu0)*Hy(i,j)+(1-exp(-sigHx(i)*dt/mu0))/(dx*sigHx(i))*(Ez(i+1,j)-Ez(i,j));
%             Hx(i,j) = Hx(i,j)-z0/2*(Ez(i,j+1)-Ez(i,j));
%             Hy(i,j) = Hy(i,j)+z0/2*(Ez(i+1,j)-Ez(i,j));

    % Update electric field 更新电场
    % Ez
    for i = 1+bpml:Nx+bpml+1
        for j = 1+bpml:Ny+bpml+1
           Ez(i,j) = CA*Ez(i,j)+CB* ...
                      ((Hy(i,j)-Hy(i-1,j))/dx-(Hx(i,j)-Hx(i,j-1))/dy);
        end
    end
    % PML
    % 左
    for i = 2:bpml
        for j = 1+bpml:Ny+bpml+1
            Ezx(i,j) = exp(-sigEx(i)*dt/eps0)*Ezx(i,j)+(1-exp(-sigEx(i)*dt/eps0))/(dx*sigEx(i))*(Hy(i,j)-Hy(i-1,j));
            Ezy(i,j) = Ezy(i,j)-1/(2*z0)*(Hx(i,j)-Hx(i,j-1));
            Ez(i,j) = Ezx(i,j)+Ezy(i,j);
        end
    end
    % 右
    for i = Nx+bpml+2:Nx+2*bpml
        for j = 1+bpml:Ny+bpml+1
            Ezx(i,j) = exp(-sigEx(i)*dt/eps0)*Ezx(i,j)+(1-exp(-sigEx(i)*dt/eps0))/(dx*sigEx(i))*(Hy(i,j)-Hy(i-1,j));
            Ezy(i,j) = Ezy(i,j)-1/(2*z0)*(Hx(i,j)-Hx(i,j-1));
            Ez(i,j) = Ezx(i,j)+Ezy(i,j);
        end
    end
    % 下
    for i = 1+bpml:Nx+bpml+1
        for j = 2:bpml
            Ezx(i,j) = Ezx(i,j)+1/(2*z0)*(Hy(i,j)-Hy(i-1,j));
            Ezy(i,j) = exp(-sigEy(j)*dt/eps0)*Ezy(i,j)-(1-exp(-sigEy(j)*dt/eps0))/(dy*sigEy(j))*(Hx(i,j)-Hx(i,j-1));
            Ez(i,j) = Ezx(i,j)+Ezy(i,j);
        end
    end
    % 上
    for i = 1+bpml:Nx+bpml+1
        for j = Ny+bpml+2:Ny+2*bpml
            Ezx(i,j) = Ezx(i,j)+1/(2*z0)*(Hy(i,j)-Hy(i-1,j));
            Ezy(i,j) = exp(-sigEy(j)*dt/eps0)*Ezy(i,j)-(1-exp(-sigEy(j)*dt/eps0))/(dy*sigEy(j))*(Hx(i,j)-Hx(i,j-1));
            Ez(i,j) = Ezx(i,j)+Ezy(i,j);
        end
    end
    % 四个角
    % 左下
    for i = 2:bpml
        for j = 2:bpml
            Ezx(i,j) = exp(-sigEx(i)*dt/eps0)*Ezx(i,j)+(1-exp(-sigEx(i)*dt/eps0))/(dx*sigEx(i))*(Hy(i,j)-Hy(i-1,j));
            Ezy(i,j) = exp(-sigEy(j)*dt/eps0)*Ezy(i,j)-(1-exp(-sigEy(j)*dt/eps0))/(dy*sigEy(j))*(Hx(i,j)-Hx(i,j-1));
            Ez(i,j) = Ezx(i,j)+Ezy(i,j);
        end
    end
    % 左上
    for i = 2:bpml
        for j = Ny+bpml+2:Ny+2*bpml
            Ezx(i,j) = exp(-sigEx(i)*dt/eps0)*Ezx(i,j)+(1-exp(-sigEx(i)*dt/eps0))/(dx*sigEx(i))*(Hy(i,j)-Hy(i-1,j));
            Ezy(i,j) = exp(-sigEy(j)*dt/eps0)*Ezy(i,j)-(1-exp(-sigEy(j)*dt/eps0))/(dy*sigEy(j))*(Hx(i,j)-Hx(i,j-1));
            Ez(i,j) = Ezx(i,j)+Ezy(i,j);
        end
    end
    % 右上
    for i = Nx+bpml+2:Nx+2*bpml
        for j = Ny+bpml+2:Ny+2*bpml
            Ezx(i,j) = exp(-sigEx(i)*dt/eps0)*Ezx(i,j)+(1-exp(-sigEx(i)*dt/eps0))/(dx*sigEx(i))*(Hy(i,j)-Hy(i-1,j));
            Ezy(i,j) = exp(-sigEy(j)*dt/eps0)*Ezy(i,j)-(1-exp(-sigEy(j)*dt/eps0))/(dy*sigEy(j))*(Hx(i,j)-Hx(i,j-1));
            Ez(i,j) = Ezx(i,j)+Ezy(i,j);
        end
    end
    % 右下
    for i = Nx+bpml+2:Nx+2*bpml
        for j = 2:bpml
            Ezx(i,j) = exp(-sigEx(i)*dt/eps0)*Ezx(i,j)+(1-exp(-sigEx(i)*dt/eps0))/(dx*sigEx(i))*(Hy(i,j)-Hy(i-1,j));
            Ezy(i,j) = exp(-sigEy(j)*dt/eps0)*Ezy(i,j)-(1-exp(-sigEy(j)*dt/eps0))/(dy*sigEy(j))*(Hx(i,j)-Hx(i,j-1));
            Ez(i,j) = Ezx(i,j)+Ezy(i,j);
        end
    end
%             Ezx(i,j) = exp(-sigEx(i)*dt/eps0)*Ezx(i,j)+(1-exp(-sigEx(i)*dt/eps0))/(dx*sigEx(i))*(Hy(i,j)-Hy(i-1,j));
%             Ezy(i,j) = Ezy(i,j)-1/(2*z0)*(Hx(i,j)-Hx(i,j-1));
%             Ez(i,j) = Ezx(i,j)+Ezy(i,j);
%             Ezx(i,j) = Ezx(i,j)+1/(2*z0)*(Hy(i,j)-Hy(i-1,j));
%             Ezy(i,j) = exp(-sigEy(j)*dt/eps0)*Ezy(i,j)-(1-exp(-sigEy(j)*dt/eps0))/(dy*sigEy(j))*(Hx(i,j)-Hx(i,j-1));
%             Ez(i,j) = Ezx(i,j)+Ezy(i,j);
         
    % Set boundary conditions 设置边界条件
    
    % Visualize fields 可视化场    
    imagesc(Ez');
    shading flat;caxis([-1.0 1.0]);axis image;axis xy; 
    title(['Ez, step ',int2str(n)]);xlabel('i');ylabel('j');
    drawnow;

end

%***********************************************************************
% END TIME-STEPPING LOOP 结束循环
%***********************************************************************