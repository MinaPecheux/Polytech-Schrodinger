% general variables
x_min = -15;
x_max = 15;
delta_x = 0.1;
x0 = -1;
sigma = 0.5;
delta_t = 0.01;
nb_steps = 14;

x_values = x_min: delta_x :x_max ;
x_values = x_values'; 
n = length(x_values);

% declaration of main matrices
laplacian = -1 / delta_x^2 * (-2*eye(n) + diag(ones(n-1,1),-1) + diag(ones(n-1,1), 1) );
V0 = zeros(n);

H = laplacian + V0;

psi = exp( - (x_values - x0).^2 / (2*sigma^2) ).^2;

plot(x_values, abs(psi / norm(psi)).^2);
hold on 

for i=1:nb_steps
%     % Explicit Euler method
%      psi = psi + delta_t * (-1j) * H * psi;
%      psi = psi / norm(psi);

% Implicit Euler method
psi = inv (eye(n) + 1j * delta_t * H) * psi; 

%      % cranck nicholson
%      psi = inv(eye(n) + 1j * delta_t/2 * H) * (eye(n) - 1j * delta_t / 2 * H) * psi;  

     if mod(i,2) == 0
         plot(x_values, abs(psi).^2, 'r')
         hold on
     end
        
end
