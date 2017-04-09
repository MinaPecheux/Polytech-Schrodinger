function oneD_TimeDependant_bis
    % general variables
    % space discrete vector
    x_min = -15;
    x_max = 15;
    delta_x = 0.1;
    x_values = (x_min:delta_x:x_max)';
    n = length(x_values);
    
    % time discrete vector
    delta_t = 0.01;
    nb_steps = 8;
    T = 0.:delta_t:delta_t*nb_steps;
    
    % declaration of main matrices
    laplacian = -1 / delta_x^2 * (-2*eye(n) + diag(ones(n-1,1),-1) + diag(ones(n-1,1), 1) );
    V0 = zeros(n);

    H = laplacian + V0;
    % functions associated with the equation
    f = @(t,x)((-1j) * H * x);
    
    % initial value of wave function
    x0 = -1.;
    sigma = 0.5;
    k = 5.;
    psi0 = exp(-(x_values - x0).^2 / (2*sigma^2) ).*exp(1j * k * x_values);
    
    % initalization of psi for all the different methods
    psi_EXP = zeros(n, nb_steps);
    psi_EXP(:,1) = psi0;

    psi_IMP = zeros(n, nb_steps);
    psi_IMP(:,1) = psi0;

    psi_CN = zeros(n, nb_steps);
    psi_CN(:,1) = psi0;

    for i = 2:nb_steps
        % Explicit Euler method
        psi_EXP(:,i) = (eye(n) - 1j * delta_t * H) * psi_EXP(:,i-1);
        psi_EXP(:,i) = psi_EXP(:,i) / norm(psi_EXP(:,i));

        % Implicit Euler method
        psi_IMP(:,i) = (eye(n) + 1j * delta_t * H) \ psi_IMP(:,i-1);
        psi_IMP(:,i) = psi_IMP(:,i) / norm(psi_IMP(:,i));

        % cranck nicholson
        psi_CN(:,i) = (eye(n) + 1j * delta_t/2 * H) \ (eye(n) - 1j * delta_t / 2 * H) * psi_CN(:,i-1);
    end
    psi_RK4 = RK4(psi0, T, f);
    
    % plotting the graphs
    figure(1)
    % Euler EXP
    subplot(2,2,1)
    hold on;
    plot(x_values, abs(psi_EXP(:,1)/norm(psi_EXP(:,1))).^2, 'DisplayName', 't = 1')
    for i = 3:2:nb_steps - 1
        plot(x_values, abs(psi_EXP(:,i)).^2, 'DisplayName', ['t = ', num2str(i)])
    end
    title('Euler explicit');
    legend('show');
    % Euler IMP
    subplot(2,2,2)
    hold on;
    plot(x_values, abs(psi_IMP(:,1)/norm(psi_IMP(:,1))).^2, 'DisplayName', 't = 1')
    for i = 3:2:nb_steps - 1
        plot(x_values, abs(psi_IMP(:,i)).^2, 'DisplayName', ['t = ', num2str(i)])
    end
    title('Euler implicit');
    legend('show');
    % Cranck Nicholson
    subplot(2,2,3)
    hold on;
    for i = 1:2:nb_steps - 1
        plot(x_values, abs(psi_CN(:,i)).^2, 'DisplayName', ['t = ', num2str(i)])
    end
    title('Cranck Nicholson');
    legend('show');
    % RK4
    subplot(2,2,4)
    hold on;
    for i = 1:2:nb_steps - 1
        plot(x_values, abs(psi_RK4(:,i)).^2, 'DisplayName', ['t = ', num2str(i)])
    end
    title('Runge Kutta 4');
    legend('show');
end

% ------------------------------------------------------------------
% Function that calculates approximate y vector with RK4 method
% Input : y0 = initial value for the function
%         T = studied time interval
%         f = function of the Cauchy problem associated to the equation
% Output : y = approximate solution vector for the equation
% ------------------------------------------------------------------
function y = RK4(y0, T, f)
    N = length(T);
    y = zeros(length(y0), N);
    y(:,1) = y0;
    for i = 1:N-1
        delta_t = T(i+1) - T(i);
        t1 = T(i);
        t2 = T(i) + delta_t/2;
        t3 = T(i) + delta_t/2;
        t4 = T(i) + delta_t;
        y1 = y(:,i);
        y2 = y1 + delta_t/2 * f(t1, y1);
        y3 = y1 + delta_t/2 * f(t2, y2);
        y4 = y1 + delta_t * f(t3, y3);
        y(:,i+1) = y(:,i) + delta_t*(1/6*f(t1,y1) + 2/6*f(t2,y2) + 2/6*f(t3,y3) + 1/6*f(t4,y4));
    end
end
