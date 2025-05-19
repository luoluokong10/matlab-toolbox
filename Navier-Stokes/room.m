% RANS (Reynolds-Averaged Navier-Stokes) using Boussinesq approximation
% Finite Differential Method for 2D incompressible flow
clear; clc; close all

% Grid setup
Lx = 1; Ly = 1;
nx = 64; ny = 64;
dx = Lx / nx; dy = Ly / ny;
x = linspace(0, Lx, nx); y = linspace(0, Ly, ny);
[X, Y] = meshgrid(x, y);

% Parameters
dt = 0.005;
tmax = 5;
nt = round(tmax / dt);
nu = 1e-3;
Cmu = 0.09;

% 初值、小涡
u = zeros(ny, nx);
v = zeros(ny, nx);
radius = 0.1;
mask = (X - 0.5).^2 + (Y - 0.5).^2 < radius^2;
u(mask) = 1;
v(mask) = -1;

% 湍流粘性场
k = 1e-2 * ones(ny, nx);
eps = 1e-2 * ones(ny, nx);
nu_t = Cmu * k.^2 ./ eps; % Boussinesq
total_nu = nu + nu_t;

% For storing energy
times = zeros(nt,1);
energy = zeros(nt,1);

for it = 1:nt
    % nu_t = Cmu * k.^2 ./ eps;
    % total_nu = nu + nu_t;

    % Laplacian
    uxx = (circshift(u, [0, -1]) - 2*u + circshift(u, [0, 1])) / dx^2;
    uyy = (circshift(u, [-1, 0]) - 2*u + circshift(u, [1, 0])) / dy^2;
    vxx = (circshift(v, [0, -1]) - 2*v + circshift(v, [0, 1])) / dx^2;
    vyy = (circshift(v, [-1, 0]) - 2*v + circshift(v, [1, 0])) / dy^2;

    % Differential
    ux = (circshift(u, [0, -1]) - circshift(u, [0, 1])) / (2*dx);
    uy = (circshift(u, [-1, 0]) - circshift(u, [1, 0])) / (2*dy);
    vx = (circshift(v, [0, -1]) - circshift(v, [0, 1])) / (2*dx);
    vy = (circshift(v, [-1, 0]) - circshift(v, [1, 0])) / (2*dy);

    adv_u = u .* ux + v .* uy;
    adv_v = u .* vx + v .* vy;

    u_diff = total_nu .* (uxx + uyy);
    v_diff = total_nu .* (vxx + vyy);

    % 迭代
    u = u + dt * (-adv_u + u_diff);
    v = v + dt * (-adv_v + v_diff);

    % 时间戳、能量
    times(it) = it * dt;
    energy(it) = sum(sum(0.5 * (u.^2 + v.^2))) * dx * dy;

    % 能量可视化
    if mod(it, round(nt/5)) == 0 || it == 1%绘制5次，每次更新前置
        plot(times(1:it), energy(1:it), 'LineWidth', 2);
        xlabel('Time'); ylabel('Kinetic Energy');
        title('Energy vs Time');
        drawnow;
    end
end

% 局部放大
zp = BaseZoom();
zp.run;