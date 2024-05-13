clear 
close all
clc

% Just using the coupled conduction solution but with very small N (but can't actually use 0)
digits(64);
T1 = 1000;
T2 = 500;
X = 10000; % Need very large X for good solution, results in X^2 matrix so careful not to exceed memory limits
N_param = 0;
tau_L = [0.01,0.1,0.5,1,2,10]; 

cur_folder = matlab.desktop.editor.getActiveFilename;
cur_folder = fileparts(cur_folder); % Just want the folder path

PM_kappa = tau_L/X;

max_itr = 10^6;
tol = 10^(-10);

% Precalculate these exponential integrals for all (i-j);
% we then map values vector indices 1:(2*X-1) to values -(X-1):(X-1) by subtracting X
% i.e., exp_int_term(1) stores the value for |i-j| = 0

theta_final = zeros(X+1,length(tau_L));
error = zeros(X+1,length(tau_L));
num_itrs = zeros(length(tau_L),1);
parfor i0 = 1:length(tau_L)
    exp_int_term = zeros(2*X-1,1); % for i = 2:X, j = 1:X+1, there are 2X-1 unique values of i-j ranging from -(X-1) to (X-1)
    exp_int_term_v2 = zeros(2*X-1,1);
    for i = 1:(2*X-1)
        exp_int_term(i) = expint(3,sym(abs(i-X)*PM_kappa(i0))) - expint(3,sym(abs(i-X+1)*PM_kappa(i0)));
        exp_int_term_v2(i) = 2*expint(3,sym(abs(i-X)*PM_kappa(i0))) - (expint(3,sym(abs(i-X+1)*PM_kappa(i0))) + expint(3,sym(abs(i-X-1)*PM_kappa(i0))));
    end
    
    theta = linspace(1,T2/T1,X+1)'; % Initial guess (theta = nondimensional temperature T/T1)
    theta(1) = 1;
    theta(end) = T2/T1;
    theta4_old = theta.^4; % theta4 = theta^4


    % Build matrix (the RHS, in modest)
    M = zeros(X+1);
    M(1,1) = 1; % Ensures inverting matrix will give theta(1)^4 = 1
    M(end,end) = 1/(T2/T1)^3; % Ensures that inverting matrix will give theta(end)^4 = (T2/T1)^4
    for i = 2:X
        for j = 2:X
            M(i,j) = exp_int_term(i-j+X)-exp_int_term(i-(j+1)+X);
            M(i,j) = exp_int_term_v2(i-j+X);
        end
        M(i,1) = -exp_int_term(i-1+X);
        M(i,X+1)=exp_int_term(i-(X+1)+X);
    end
    M = M;

    LHS = zeros(X+1,1); % constant values in the LHS
    LHS(1) = 1; % LHS for i = 1 is always T1/T1 = 1
    LHS(end) = T2/T1; %
    if N_param > 0
        under_relax = (0.1)/(max(N_param,10^(-6))*X^2)*tau_L(i0); % Through guess and checking found that this works well
    
        for i = 2:X
            LHS(i) = (theta(i-1)-2*theta(i)+theta(i+1))*(2*N_param)/PM_kappa(i0);
        end

        for n = 1:max_itr
            
            theta4 = M\LHS; % Solve M*X = LHS
            theta4 = theta4_old + under_relax*(theta4-theta4_old); % Use underrelaxation to get updated theta4
            theta4_old = theta4; % Update theta4_old
            theta_old = theta; % Update theta_old
            theta = theta4.^(1/4); % Update theta

            for i = 2:X
                LHS(i) = (theta(i-1)-2*theta(i)+theta(i+1))*(2*N_param)/PM_kappa(i0);
            end
            
            error(:,i0) = M*theta4 - LHS;
            
            if max(abs(error(:,i0))) < tol
                num_itrs(i0) = n;
                break;
            end
        
        end
    else
        theta4 = M\LHS;
        theta = theta4.^(1/4);
    end
    theta_final(:,i0) = theta;
end

figure
hold on;

legend_str = cell(1,length(tau_L));
T = theta_final*T1;

nondim_power = ((theta_final*T1).^4-T2^4)/(T1^4-T2^4); % Nondimensional emissive power

for i0 = 1:length(tau_L)
    plot(linspace(0,1,X+1),nondim_power(:,i0));
    legend_str{i0} = char(strcat("\tau_L = ",string(tau_L(i0))));
end
legend(legend_str)
legend_str = [{'Non-dim position'},legend_str];
varTypes = cell(1,length(tau_L));
for i0 = 1:(1+length(tau_L))
    varTypes{i0} = 'double';
end
results_table = table('Size',[X+1,1+length(tau_L)],'VariableTypes',varTypes,'VariableNames',legend_str);

results_table{:,1} = linspace(0,1,X+1)';
results_table{:,2:(length(tau_L)+1)} = nondim_power;

digits(32);
save(strcat(cur_folder,"\ParallelPlatesPM_Exact"),'results_table')