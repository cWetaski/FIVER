%% Sphere_SA
% Test the accuracy of the surface normal and surface area estimation for a voxelized sphere.
clear
close all
clc
if isempty(gcp('nocreate'))
    parpool("Threads",8);
end
%% Params
format shortG
D = 10:100; % vx
ns = 1:30;
N_tests = length(D);
N_ns = length(ns);
plot_results = true;

%% Analytic Solution
SA = 4*pi.*(D./2).^2; % [vx^2]: Surface area of sphere
V = 4/3*pi.*(D./2).^3; % [vx^3]: Volume of sphere


SA_est = zeros(N_tests,1);
SA_error = zeros(N_tests,1);
V_est = zeros(1,N_tests);
V_error = zeros(1,N_tests);
VS_surf_norms_saved = [];
parfor i = 1:N_tests
    disp(D(i))
    %% Generate Voxel Space
    d = D(i);
    VS_opaq = generateSphere(d);
    % Evaluate error
    V_est(i) = sum(VS_opaq,'all');
    V_error(i) = (V_est(i)-V(i))/V(i)*100;
    VS_opaq_outer = ~VS_opaq;
    
    [VS_surf_norms, VS_surf_areas, ~] = getNormalsAndSurfaceAreas(VS_opaq,ns(1));
    [VS_surf_norms_outer,VS_surf_areas_outer,~] = getNormalsAndSurfaceAreas(VS_opaq_outer,ns(1));
    SA_est_outer(i) = sum(VS_surf_areas_outer,'all');
    SA_error_outer(i) = (SA_est_outer(i)-SA(i))/SA(i)*100;
    SA_est(i) = sum(VS_surf_areas,"all");
    SA_error(i) = (SA_est(i)-SA(i))/SA(i)*100;
    SA_avg(i) = (SA_est(i)+SA_est_outer(i))/2;
    SA_error_avg(i) = (SA_avg(i)-SA(i))/SA(i)*100;
end
if plot_results
    cmap = load("PREC.mat").color_scheme;
    f = figure;
    hold on
    p = plot(D,SA_error,'Color',cmap(1),'LineWidth',1);
    p_V = plot(D,V_error,'Color','black','LineWidth',1);
    legend(gca,[p,p_V],{"Surface area","Volume"},'Interpreter','latex')
    xlabel('Diameter (vx)','Interpreter','latex');
    ylabel('Error (\%)','Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');
end
disp("FINISHED DIAMETER NOW DOING NS")
VS_opaq = generateSphere(D(end));
cen = D(end)/2;
theta_error = zeros(N_ns,1);
count = 0;
X = [];
Y = [];
Z = [];
U = [];
V = [];
W = [];
d = D(end);
SA_est_ns = zeros(N_ns,1);
SA_error_ns = zeros(N_ns,1);
for j = 1:N_ns
    disp(ns(j));
    [VS_surf_norms, VS_surf_areas, ~] = getNormalsAndSurfaceAreas(VS_opaq,ns(j));
    for xi = 1:d+2
        for yi = 1:d+2
            for zi = 1:d+2
                surf_norm = VS_surf_norms{xi,yi,zi};
                if ~isempty(surf_norm)
                    count = count + 1;
                    X(count) = xi-0.5;
                    Y(count) = yi-0.5;
                    Z(count) = zi-0.5;
                    U(count) = surf_norm(1);
                    V(count) = surf_norm(2);
                    W(count) = surf_norm(3);
                end
            end
        end
    end
    P = [X',Y',Z'];
    exact_lines = P-cen;
    U_lines = exact_lines(:,1);
    V_lines = exact_lines(:,2);
    W_lines = exact_lines(:,3);
    X_lines = cen*ones(count,1);
    Y_lines = cen*ones(count,1);
    Z_lines = cen*ones(count,1);
    SA_est_ns(j) = sum(VS_surf_areas,'all');
    SA_error_ns(j) = (SA_est_ns(j)-SA(end))/SA(end)*100;
    exact_norms = exact_lines./vecnorm(exact_lines,2,2);
    theta_error(j) = mean(real(acos(dot(exact_norms,[U',V',W'],2)))*180/pi,'all');
end
figure
plot(ns,SA_error_ns,'Color',cmap(1),'LineWidth',1);
xlabel('Surface Normal Neighborhood Size (vx)','Interpreter','latex');
ylabel('Surface Area Error (\%)','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex');
ylim([floor(min(SA_error_ns)),max(0,ceil(max(SA_error_ns)))]);

disp(SA_err2)
fprintf("Hollow Sphere Surface Area \n");
fprintf("Analytic: %0.3f \n",SA(end,1));
fprintf("Voxelized Estimation: %0.3f \n", SA_est(end,1));
fprintf("Error: %0.3f %% \n",SA_error(end,1));

figure
yyaxis left
plot(ns,theta_error)
xlabel('Surface Normal Neighborhood Size (vx)','Interpreter','latex');
ylabel('Mean Angular Error $\bar{\epsilon}_\theta$ ($^{\circ}$)','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex');
ylim([0,ceil(max(theta_error))]);
yyaxis right
plot(ns,SA_error_ns,'Color',cmap(2),'LineWidth',1);
ylabel('Surface Area Error (\%)','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex');
ylim([floor(min(SA_error_ns)),max(0,ceil(max(SA_error_ns)))]);

function VS_opaq = generateSphere(d)
    VS_opaq = false(d+2,d+2,d+2);
    r = d/2;
    cen = (d+2)/2;
    for xi = 1:(d+2)
        x_cen = xi-0.5;
        rx = (cen-x_cen)^2;
        for yi = 1:(d+2)
            y_cen = yi-0.5;
            ry = (cen-y_cen)^2;
            for zi = 1:(d+2)
                z_cen = zi-0.5;
                rz = (cen-z_cen)^2;
                if rx+ry+rz < r*r
                    VS_opaq(xi,yi,zi) = true;
                end
            end
        end
    end
end