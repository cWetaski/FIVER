clear 
close all
clc

T1 = 800;
T2 = 400;
X = 1600;
tau_L = 1;

N_param = [10,1,0.1,0.01,0.001];



PM_kappa = tau_L/X;

max_itr = 10^5;
tol = 10^(-8);

theta_i = zeros(X+1,length(N_param));
nondim_x = linspace(0,1,X+1)';

% Get current folder
cur_folder = matlab.desktop.editor.getActiveFilename;
cur_folder = fileparts(cur_folder); % Just want the folder path

try % Load existing results (this way we can improve on existing values instead of starting from scratch)
    results_table = load(strcat(cur_folder,"ParallelPlatesPM_Cond_Exact.mat")).results_table;
    x_table = results_table{:,1};
    for i = 1:length(N_param)
        string_N = strcat("N = ",string(N_param(i)));
        try 
            theta_table = results_table{:,string_N};
            if length(x_table) == length(nondim_x)
                theta_i(:,i) = theta_table;
            else
                theta_i(:,i) = interp1(x_table,theta_table,nondim_x);
            end
        catch
            theta_i(:,i) = linspace(1,T2/T1,X+1)';
        end
    end
catch e
    for i = 1:length(N_param)
        theta_i(:,i) = linspace(1,T2/T1,X+1)';
    end
end


% Precalculate these exponential integrals for all (i-j);
exp_int_term = zeros(X,1); % for i = 2:X, j = 1:X+1, there are 2X-1 unique values of i-j ranging from -(X-1) to (X-1)
% we then map values vector indices 1:(2*X-1) to values -(X-1):(X-1) by subtracting X
% i.e., exp_int_term(1) stores the value for |i-j| = 0
for i = 1:(2*X-1)
exp_int_term(i) = double(expint(3,sym(abs(i-X)*PM_kappa)) - expint(3,sym(abs(i-X+1)*PM_kappa)));
end


theta_final = zeros(X+1,length(N_param));
error_val = zeros(length(N_param),1);
num_itrs = zeros(length(N_param),1);

% Build LHS matrix
LHS_M = spalloc(X+1,X+1,(X+1)*3-4);
LHS_M = LHS_M - speye(X+1)*2; % Set diagonals to -2
for i = 2:X 
    for j = (i-1):2:(i+1) % Only 2 values for row
        LHS_M(i,j) = 1;
    end
end
LHS_M(1,1) = 1; % These values are just equal to to the RHS
LHS_M(end,end) = 1;


% Build RHS matrix
RHS_M = zeros(X+1);

for i = 2:X
    for j = 2:X
        RHS_M(i,j) = exp_int_term(i-j+X)-exp_int_term(i-(j+1)+X);
    end
    RHS_M(i,1) = -exp_int_term(i-1+X);
    RHS_M(i,X+1)=exp_int_term(i-(X+1)+X);
end
RHS_M = RHS_M*PM_kappa/(2*N_param(1));

for i0 = 1:length(N_param)
    
    if i0 > 1
        RHS_M = RHS_M*N_param(i0-1)/N_param(i0); % update RHS
    end
    RHS_M(1,1) = 1; % Ensures inverting matrix will give theta(1)^4 = 1
    RHS_M(end,end) = 1/(T2/T1)^3; % Ensures that inverting matrix will give theta(end)^4 = (T2/T1)^4

    if N_param(i0) >= 0.001 % Conduction dominates
        under_relax = min(0.5,0.5*N_param(i0)*10); % under relaxation factor for conduction dominating
        theta = theta_i(:,i0); % Initial guess
        theta_old = theta;
    
        RHS = RHS_M*theta.^4;
        error_val(i0) = max(abs(LHS_M*theta-RHS)); % intial error
        n = 0;
        while error_val(i0) > tol && n < max_itr        
            n = n+1;

            theta = LHS_M\RHS; % Solve LHS_M*X = RHS
            theta = theta_old + under_relax*(theta-theta_old);
            theta_old = theta;

            RHS = RHS_M*theta.^4; % Update RHS

            error_val(i0) = max(abs(LHS_M*theta-RHS));
            
            if mod(n,1) == 0
                plot(linspace(0,1,X+1),theta);drawnow();
                disp(error_val(i0));
            end
        end
    else % Radiation dominates (sufficiently under relaxed solution in the other direction is usually better regardless)
        under_relax = (0.1)/(N_param(i0)*X^2); % Through guess and checking found that this works well
        theta = theta_i(:,i0);
        theta(1) = 1;
        theta(end) = T2/T1;
        theta4_old = theta.^4;
        theta4 = theta4_old;

        LHS = LHS_M*theta;
        error_val(i0) = max(abs(RHS_M*theta4-LHS)); % initial error;
        n = 0;
        while error_val(i0) > tol && n < max_itr
            n = n+1;

            theta4 = RHS_M\LHS; % Solve RHS_M*X = LHS
            theta4 = theta4_old + under_relax*(theta4-theta4_old); % Use underrelaxation to get updated theta4
            theta4_old = theta4; % Update theta4_old
            theta = theta4.^(1/4); % Update theta
            
            LHS = LHS_M*theta; % Update LHS

            error_val(i0) = max(abs(RHS_M*theta4 - LHS)); % Check error
            if mod(n,1000) == 0
                plot(linspace(0,1,X+1),theta);drawnow();
                disp(error_val(i0))
            end
        end
    end
    num_itrs(i0) = n;
    theta_final(:,i0) = theta;
    
end

figure
hold on;

legend_str = cell(1,length(N_param));
for i = 1:length(N_param)
    plot(linspace(0,1,X+1),theta_final(:,i));
    legend_str{i} = char(strcat("N = ",string(N_param(i))));
end
legend(legend_str)
legend_str = [{'Non-dim position'},legend_str];
varTypes = cell(1,length(N_param));
for i = 1:(1+length(N_param))
    varTypes{i} = 'double';
end
results_table = table('Size',[X+1,1+length(N_param)],'VariableTypes',varTypes,'VariableNames',legend_str);

results_table{:,1} = linspace(0,1,X+1)';
results_table{:,2:(length(N_param)+1)} = theta_final;

save(strcat(cur_folder,"\ParallelPlatesCondPM_Exact"),'results_table');