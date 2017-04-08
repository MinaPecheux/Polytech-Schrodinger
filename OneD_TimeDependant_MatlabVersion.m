function OneD_TimeDependant_MatlabVersion()

% general variables
x_min = -15;
x_max = 15;
delta_x = 0.1;
x0 = -1;
sigma = 0.5;
delta_t = 0.01;
nb_steps =18;

x_values = (x_min: delta_x :x_max)' ;
n = length(x_values);

% declaration of main matrices
laplacian = -1 / delta_x^2 * (-2*eye(n) + diag(ones(n-1,1),-1) + diag(ones(n-1,1), 1) );
V0 = zeros(n);

H = laplacian + V0;

k = 0;
psi = exp( - (x_values - x0).^2 / (2*sigma^2) ).*exp(1j * k * x_values);

% initalization of psi for all the different methods
psi_EXP= zeros(n, nb_steps);
psi_EXP(:,1) = psi;

psi_IMP = zeros(n, nb_steps);
psi_IMP(:,1) = psi;

psi_CN = zeros(n, nb_steps);
psi_CN(:,1) = psi;

psi_RK4 = zeros(n, nb_steps);
psi_RK4(:,1) = psi;


for i=2:nb_steps
    
    % Explicit Euler method
    psi_EXP(:,i) = psi_EXP(:,i-1) + delta_t * (-1j) * H * psi_EXP(:,i-1); 
    psi_EXP(:,i) = psi_EXP(:,i) / norm(psi_EXP(:,i));
    
    % Implicit Euler method
    psi_IMP(:,i) = (eye(n) + 1j * delta_t * H) \ psi_IMP(:,i-1);
    
    % cranck nicholson
    psi_CN(:,i) = (eye(n) + 1j * delta_t/2 * H) \ (eye(n) - 1j * delta_t / 2 * H) * psi_CN(:,i-1);
    
    %RK4
    k1 = f( H, psi_RK4(:,i-1)                      );
    k2 = f( H, psi_RK4(:,i-1) + (delta_t / 2) * k1 );
    k3 = f( H, psi_RK4(:,i-1) + (delta_t / 2) * k2 );
    k4 = f( H, psi_RK4(:,i-1) +  delta_t      * k3 );
    
    psi_RK4(:,i) = psi_RK4(:,i-1) + (delta_t / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
   % psi_RK4 = psi_RK4 / norm(psi_RK4);  
        
end

%plotting the graphs
for i=1:nb_steps
    
    if mod(i,2) == 0
        subplot(2,2,1)
        plot(x_values, abs( psi_EXP(:,i) ).^2, 'r')
        legend('Explicite')
        hold on
        
        subplot(2,2,2)
        plot(x_values, abs( psi_IMP(:,i) ).^2, 'k')
        legend('Implicite')
        hold on
        
        subplot(2,2,3)
        plot(x_values, abs( psi_CN(:,i) ).^2, 'b')
        legend('Cranck-Nicholson')
        hold on
        
        subplot(2,2,4)
        plot(x_values, abs( psi_RK4(:,i) ).^2, 'm')
        legend('RK4')
        hold on
        
    end
end

end

function f = f(H, psi)

f = (-1j) * H * psi;

end

