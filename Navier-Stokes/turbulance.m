function turbulance
    % 网格和时间设置
    nx = 64; ny = 64;
    Lx = 1; Ly = 1;
    dx = Lx / nx; dy = Ly / ny;
    dt = 0.005; tf = 2;
    nt = round(tf / dt);
    
    % 网格生成
    [x, y] = meshgrid(linspace(0, Lx, nx), linspace(0, Ly, ny));
    
    % 初始速度场：随机扰动
    u = 0.1 * (2*rand(ny, nx) - 1); % x方向
    v = 0.1 * (2*rand(ny, nx) - 1); % y方向
    
    % 压力场初始化
    p = zeros(ny, nx);

    % 模拟两个不同的雷诺数
    Re_vals = [100, 1000]; % 一个层流，一个趋近湍流
    figure;
    
    for i = 1:length(Re_vals)
        u1 = u; v1 = v; p1 = p;
        Re = Re_vals(i);
        
        for t = 1:nt
            % 中心差分梯度和拉普拉斯
            uxx = (circshift(u1, [0 -1]) - 2*u1 + circshift(u1, [0 1])) / dx^2;
            uyy = (circshift(u1, [-1 0]) - 2*u1 + circshift(u1, [1 0])) / dy^2;
            vxx = (circshift(v1, [0 -1]) - 2*v1 + circshift(v1, [0 1])) / dx^2;
            vyy = (circshift(v1, [-1 0]) - 2*v1 + circshift(v1, [1 0])) / dy^2;

            % 非线性对流项
            u_adv = u1 .* (circshift(u1, [0 -1]) - circshift(u1, [0 1])) / (2*dx) + ...
                    v1 .* (circshift(u1, [-1 0]) - circshift(u1, [1 0])) / (2*dy);
            v_adv = u1 .* (circshift(v1, [0 -1]) - circshift(v1, [0 1])) / (2*dx) + ...
                    v1 .* (circshift(v1, [-1 0]) - circshift(v1, [1 0])) / (2*dy);

            % 更新速度（显式）
            u1 = u1 + dt * (-u_adv + (1/Re) * (uxx + uyy));
            v1 = v1 + dt * (-v_adv + (1/Re) * (vxx + vyy));
        end

        % 绘图：速度大小和方向
        subplot(1,2,i)
        speed = sqrt(u1.^2 + v1.^2);
        quiver(x, y, u1, v1, 'k')
        hold on
        contourf(x, y, speed, 20, 'LineStyle','none')
        colorbar
        title(sprintf('Re = %d', Re))
        axis equal tight
        xlabel('x'), ylabel('y')
    end
end
