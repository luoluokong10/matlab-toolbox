% RANS (Reynolds-Averaged Navier-Stokes) using Boussinesq approximation
% Finite Volume Method for 2D incompressible flow
clear; clc;

% Grid setup
Lx = 1; Ly = 1;
nx = 64; ny = 64;
dx = Lx / nx; dy = Ly / ny;
x = linspace(0, Lx, nx); y = linspace(0, Ly, ny);
[X, Y] = meshgrid(x, y);

% Time parameters
dt = 0.005;
tmax = 5;
nt = round(tmax / dt);

% Physical parameters
nu = 1e-3;    % Molecular viscosity
Cmu = 0.09;   % Empirical constant in Boussinesq model

% Initial condition: small circular vortex
u = zeros(ny, nx);
v = zeros(ny, nx);
radius = 0.1;
mask = (X - 0.5).^2 + (Y - 0.5).^2 < radius^2;
u(mask) = 1;
v(mask) = -1;

% Eddy viscosity field (turbulent viscosity)
k = 1e-2 * ones(ny, nx);   % Turbulent kinetic energy
eps = 1e-2 * ones(ny, nx); % Dissipation rate
nu_t = Cmu * k.^2 ./ eps; % Boussinesq model
total_nu = nu + nu_t;

% For storing energy
times = zeros(nt,1);
energy = zeros(nt,1);

for it = 1:nt
    % Effective viscosity recomputed at each step (optional update)
    nu_t = Cmu * k.^2 ./ eps;
    total_nu = nu + nu_t;

    % Central difference for Laplacian
    uxx = (circshift(u, [0, -1]) - 2*u + circshift(u, [0, 1])) / dx^2;
    uyy = (circshift(u, [-1, 0]) - 2*u + circshift(u, [1, 0])) / dy^2;
    vxx = (circshift(v, [0, -1]) - 2*v + circshift(v, [0, 1])) / dx^2;
    vyy = (circshift(v, [-1, 0]) - 2*v + circshift(v, [1, 0])) / dy^2;

    % First-order upwind advection (simple for demonstration)
    ux = (circshift(u, [0, -1]) - circshift(u, [0, 1])) / (2*dx);
    uy = (circshift(u, [-1, 0]) - circshift(u, [1, 0])) / (2*dy);
    vx = (circshift(v, [0, -1]) - circshift(v, [0, 1])) / (2*dx);
    vy = (circshift(v, [-1, 0]) - circshift(v, [1, 0])) / (2*dy);

    adv_u = u .* ux + v .* uy;
    adv_v = u .* vx + v .* vy;

    % Diffusion term (Boussinesq with total viscosity)
    u_diff = total_nu .* (uxx + uyy);
    v_diff = total_nu .* (vxx + vyy);

    % Time integration (Euler forward)
    u = u + dt * (-adv_u + u_diff);
    v = v + dt * (-adv_v + v_diff);

    % Store energy
    times(it) = it * dt;
    energy(it) = sum(sum(0.5 * (u.^2 + v.^2))) * dx * dy;

    % Visualization
    if mod(it, round(nt/5)) == 0 || it == 1
        subplot(1,2,1);
        imagesc(x, y, sqrt(u.^2 + v.^2)); axis image; colorbar;
        title(['Velocity magnitude, t = ' num2str(it*dt)]);
        xlabel('x'); ylabel('y');

        subplot(1,2,2);
        plot(times(1:it), energy(1:it), 'LineWidth', 2);
        xlabel('Time'); ylabel('Kinetic Energy');
        title('Energy vs Time');
        drawnow;
    end
end