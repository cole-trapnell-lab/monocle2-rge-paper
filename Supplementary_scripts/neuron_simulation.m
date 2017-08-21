%branching of neural development
%the regulatory network is mannually collected
%use Huang Sui's SDEs approaches
%Author: Xiaojie Qiu Email: Tuexy.xiaojie@gmail.com 
%version 1.0

% def = struct(...
%         'a_s', 2.2, ...
%         'a_e', 6, ...
%         'a', 4, ...
%         'k',1, ...
%         'namada', 0.25, ...
%         'namada_m', 0.125, ...
%         'n',4, ...
%         'm',0.01, ...
%         'N', 14);
% 
% %check the input and set the parameters of the SDEs
% if ~nargin 
%     Parameters = def;
% elseif ~isstruct(Parameters);
%     error('MATLAB:hierarcy_attractor', 'You need to input a struct');
% else 
%     fs = {'a_s', 'a_e', 'a', 'k', 'namada', 'namada_m', 'n', ...
%         'm', 'N'};
%     for n = 1:length(fs)
%         if ~isfield(Parameters, fs{n}), Parameters.(fs{n}) = def.(fs{nm}); end
%     end
% end
% 
% %main setting
% a_s = Parameters.a_s;
% a_e = Parameters.a_e;
% a = Parameters.a;
% k = Parameters.k;
% eta = Parameters.namada;
% eta_m = Parameters.namada_m;
% n = Parameters.n; 
% m = Parameters.m;
% N = Parameters.N;

clear all;
mature_mu = 0; %0.40; 
n = 4; N = 13; a = 4; a_s = 2.2; a_e = 6; k = 1; eta = 0.25; eta_m = 0.125; eta_b = 0.1;
mx = 10;
%n = 8; N = 12; a = 8; a_s = 4.4; a_e = 12; k = 2; eta = 0.5; eta_m = 0.25; eta_b = 0.2; 
N_end = 100;
steps = 2000;
dt = N_end / steps;
cells = 200;     
mu = zeros(steps, N -1); %1000steps,N dimentions
syms x_1 x_2 x_3 x_5 x_6 x_7 x_8 x_9 x_10 x_11 x_12 x_14 x_15; 
%Aldh1L: no repression to stat3;  %%%new: clear the auto of mash1 and Hes5
%Tuj1: a_e; 
%Pax6
%F = [a_s * (x_1^n) / (1 + x_1^n + x_9^n + x_6 ^n + x_12 ^n ) - k * x_1 ; ... %add noise
F = [a_s * 1 / (1 + eta^n * (x_6 + x_12 + x_9)^n * x_15 ^n) - k * x_1 ; ... %add noise
%Mash1
%a * (eta ^ n * x_1 ^ n + x_2 ^ n + eta_b ^ n * x_11 ^ n) / (1 + eta ^ n * x_1 ^ n + x_2^n + x_7^n + eta_b ^ n * x_11 ^ n) - k * x_2 ; ...
a * (x_1 ^ n) / (1 + x_1 ^ n + x_7^n) - k * x_2 ; ...
%Brn2
a * (x_2^n) / (1 + x_2^n) - k * x_3 ; ... %no Olig2 inhibition
%a * (eta ^ n * x_2^n) / (1 + eta ^ n * x_2^n + x_9^n) - k * x_3 ; ...
%%NeuroD
%a * ((x_2^n) + x_4^n) / (1 + (x_2^n) + x_4^n) - k * x_4 ; ...
%Zic1
a * (x_2^n) / (1 + x_2^n)  - k * x_5 ; ...
%a * (eta ^ n * x_2^n) / (1 + eta ^ n * x_2^n)  - k * x_5 ; ...
%Tuj1
a_e * (x_3^n + x_5^n + x_11^n) / (1 + x_3^n + x_5^n + x_11^n) - k * x_6 ; ...
%Hes5
a * (x_1^n) / (1 +  x_1^n + x_2^n)  - k * x_7 ; ...
%Scl
a_e * (eta ^ n * x_7^n) / (1 + eta ^ n * x_7^n + x_9^n) - k * x_8; ...
%Olig2
a_e * ((eta * x_7) ^n ) / (1 + (x_8)^n + (eta * x_7) ^n  ) - k * x_9; ... %no auto no brn2 inhibition
%Stat3
a * (eta^n * x_7^n * x_8^n) / (1 + eta^n *  x_7^n * x_8^n ) - k * x_10; ...
%a * (eta_m ^ n * x_7^n * eta ^ n * x_8^n) / (1 + eta_m ^ n * x_7^n * eta ^ n * x_8^n ) - k * x_10; ...
%Myt1L
a * (x_9^n) / (1 + x_9^n) - k * x_11; ...       
%a * (eta ^ n * x_9^n) / (1 + eta ^ n * x_9^n) - k * x_11; ...       
%Aldh1L
a_e * (x_10^n) / (1 + x_10^n) - k * x_12; ...
%%Dnmt1
%a * 1 / (1 + x_10^n) - k * x_13; ...
%Sox8
a * (eta_m ^n*  x_9^n) / (1 + eta_m ^n  *x_9^n) - k * x_14; ...
mature_mu * (1 - x_15 / mx);
];

res = solve(F,x_1,x_2,x_3,x_5,x_6,x_7,x_8,x_9,x_10,x_11, x_12, x_14, x_15);

%EM method
%cell_simulate = zeros(50, 11, 1000); %1000 cells * 11 genes * 5000 steps
cell_simulate = zeros(N, steps + 1, cells);
f_em = zeros(N, steps+1);
f_init = zeros(N, 1);
f_init(1) = 2; %large expression level of Pdx6;
%f_temp = f_init;
D(N -1, N -1, steps)=0;
f_em(:,1) = f_init;
for cell = 1:cells
    f_temp = f_init;
    d = randn(1, N -1, steps) / 10;
    for sh = 1:steps
        D(:,:,sh) = (d(:,:,steps))' * d(:,:,steps);   %diffusion matrix
    end
    noise = mvnrnd(mu, 2*D, steps); %use mvnrnd to set the noise and multiply the sqrt of time
    for i = 1:steps
        f_temp = f_temp + double(subs(F, {x_1, x_2, x_3,x_5, x_6, x_7, x_8, x_9, x_10,x_11, x_12, x_14 x_15}, f_temp'))*dt + [noise(i, :) 0]' * sqrt(dt); 
        f_temp(f_temp <0) = 0; %keep non-negative of the expression value
        f_em(:, (i+1)) = f_temp;
    end
   % cell_simulate(cells, :, :) = f_em;
   cell_simulate(:, :, cell) = f_em;
end
save cell_simulate.mat; %save the data for future use;

%expression(steps+1)=0;
%figure 3 in the paper
 figure(2);
subplot(2, 2, 1);  %1st branching
x = 0:dt:N_end;
for cell = 1:cells
% plot(x, cell_simulate(1, :, cell), x, cell_simulate(2,:,cell), x, cell_simulate(6, :, cell),... %2, 3, 10????
% 'LineWidth', 2, 'Color', 'r', 'Color', 'm', 'Color', 'b')
plot(x, cell_simulate(1, :, cell), ... %2, 3, 10????
'LineWidth', 2, 'Color', 'r')
hold on;
plot(x, cell_simulate(2, :, cell), ... %2, 3, 10????
'LineWidth', 2, 'Color', 'm')
hold on;
plot(x, cell_simulate(6, :, cell),... %2, 3, 10????
'LineWidth', 2, 'Color', 'b')
hold on;
end
axis([0 N_end   0 5]);
legend('Pax6', 'Mash1', 'Hes5') %here?
ylabel('Gene expression');
grid on;
title('lst branching, Mash1 - Hes5');
subplot(2, 2, 3); %2st branching
for cell=1:cells
plot(x, cell_simulate(5,:,cell), 'y', x, cell_simulate(7,:,cell), 'c', x, cell_simulate(8,:,cell), 'g', x, cell_simulate(11,:,cell), 'b',x, cell_simulate(12,:,cell), 'm', 'LineWidth', 2);
hold on;
end
axis([0 N_end  0 7]);
title('2nd branching, Scl - Olig2');
legend('Tuj1', 'Scl', 'Olig2', 'Aldh1L', 'Sox8');
ylabel('Gene expression'); grid on;
subplot(2, 2, 2);
for cell=1:cells
plot(cell_simulate(2,1,cell), cell_simulate(6,1,cell),'--rs', 'MarkerEdgeColor','k', 'MarkerFaceColor','g', 'MarkerSize',15); %start point
plot(cell_simulate(2,:,cell), cell_simulate(6,:,cell), 'r', 'LineWidth', 2);
plot(cell_simulate(2,end,cell), cell_simulate(6,end,cell), '--rs', 'MarkerEdgeColor','k', 'MarkerFaceColor','g', 'MarkerSize',15); %start point
hold on;
end
xlabel('Mash1'), ylabel('Hes5'), grid on; axis([0 5 0 5]);
subplot(2, 2, 4);
for cell=1:cells
plot(cell_simulate(7,1,cell), cell_simulate(8,1,cell), '--rs', 'MarkerEdgeColor','k', 'MarkerFaceColor','g', 'MarkerSize',15); %start point
plot(cell_simulate(7,:,cell), cell_simulate(8,:,cell), 'r', 'LineWidth', 2);
plot(cell_simulate(7,end,cell), cell_simulate(8,end,cell), '--rs', 'MarkerEdgeColor','k', 'MarkerFaceColor','g', 'MarkerSize',15); %end point
hold on;
end
xlabel('Scl'), ylabel('Olig2'), grid on; axis([0 5 0 5]);
