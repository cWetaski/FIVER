clear
close all
clc;
c1_switch = false;
v = [0.6995,0.6784,0.361,0.5493,0.312,53.4];
data = readtable('square_laminar_Nu_data.csv');

x_vals = data.x(2:end);
Nu_x = data.Nu_x(2:end);
Nu_m = data.Nu_m(2:end);
N_vals = length(x_vals);
Nu_x_fun = zeros(N_vals,1);
Nu_m_fun = zeros(N_vals,1);
for i = 1:N_vals
    Nu_x_fun(i) = getNu_x_lam(x_vals(i),v);
    Nu_m_fun(i) = getNu_m_lam(0,x_vals(i),v);
end

v2 = fminsearch(@getMSE,v);

Nu_x_fun2 = zeros(N_vals,1);
Nu_m_fun2 = zeros(N_vals,1);
for i = 1:N_vals
    Nu_x_fun2(i) = getNu_x_lam(x_vals(i),v2);
    Nu_m_fun2(i) = getNu_m_lam(0,x_vals(i),v2);
end

[error_Nu_x,error_Nu_m] = getErrors(v);
[error_Nu_x2,error_Nu_m2] = getErrors(v2);
errors1 = [error_Nu_x;error_Nu_m];
errors2 = [error_Nu_x2;error_Nu_m2];

max(errors1)
max(errors2)
N_full = 1000;
x_full = 10.^(linspace(-3,0,N_full));

Nu_x_full = zeros(N_full,1);
Nu_m_full = zeros(N_full,1);
Nu_x_full2 = zeros(N_full,1);
Nu_m_full2 = zeros(N_full,1);

for i = 1:N_full
    Nu_x_full(i) = getNu_x_lam(x_full(i),v);
    Nu_m_full(i) = getNu_m_lam(0,x_full(i),v);
    Nu_x_full2(i) = getNu_x_lam(x_full(i),v2);
    Nu_m_full2(i) = getNu_m_lam(0,x_full(i),v2);
end

figure()
hold on
scatter(x_vals,Nu_x,'k','Marker','o');
scatter(x_vals,Nu_m,'k','Marker','o');
scatter(x_vals,Nu_x_fun,'magenta','Marker','+')
scatter(x_vals,Nu_m_fun,'magenta','Marker','+')
set(gca,'YScale','log')
set(gca,'XScale','log')


scatter(x_vals,Nu_x_fun2,'blue','Marker','x')
scatter(x_vals,Nu_m_fun2,'blue','Marker','x')
plot(x_full,Nu_x_full,'--','Color','k')
plot(x_full,Nu_m_full,'--','Color','k')
plot(x_full,Nu_x_full2,':','Color','k')
plot(x_full,Nu_m_full2,':','Color','k')

function [error_Nu_x,error_Nu_m] = getErrors(v)
    data = readtable('square_laminar_Nu_data.csv');
    x_vals = data.x(2:end);
    Nu_x = data.Nu_x(2:end);
    Nu_m = data.Nu_m(2:end);
    N_vals = length(x_vals);
    Nu_x_fun = zeros(N_vals,1);
    Nu_m_fun = zeros(N_vals,1);
    for i = 1:N_vals
        Nu_x_fun(i) = getNu_x_lam(x_vals(i),v);
        Nu_m_fun(i) = getNu_m_lam(0,x_vals(i),v);
    end
error_Nu_x = abs(Nu_x-Nu_x_fun)./Nu_x;
error_Nu_m = abs(Nu_m-Nu_m_fun)./Nu_m;
end
function MSE = getMSE(v)
[error_Nu_x,error_Nu_m] = getErrors(v);
MSE = mean(([error_Nu_x;error_Nu_m]).^2);
end

function Nu_x = getNu_x_lam(x_star,v)
    % Correlation fit to data from Table 3.15, Bhatti and Shah 1987
    a1 = 0.6365; b1 = 0.8880; c1 = 0.3702;
    a2 = 1.181; b2 = 0.1260; c2 = 56.77;
    a1 = v(1); b1 = v(2); c1 = v(3);
    a2 = v(4); b2 = v(5); c2 = v(6);
    if x_star <= 1/80 % 
        Nu_x = a1*x_star.^(-c1)+b1;
    else
        Nu_x = 2.975+a2.*x_star.^(-b2).*exp(-c2*x_star);
    end
end
function Nu_m = getNu_m_lam(x_star_L,x_star_U,v)
    % x_star = 1/Gz
    % x_star_L < x_star_U
    a1 = 0.6365; b1 = 0.8880; c1 = 0.3702;
    a2 = 1.181; b2 = 0.1260; c2 = 56.77;
    a1 = v(1); b1 = v(2); c1 = v(3);
    a2 = v(4); b2 = v(5); c2 = v(6);
    if x_star_U <= 1/80
        Nu_m = integral(@(x) a1*x.^(-c1)+b1,x_star_L,x_star_U);
    elseif x_star_L >= 1/80
        Nu_m = integral(@(x) 2.975+a2.*x.^(-b2).*exp(-c2*x),1/1000,x_star_U);
    else % x_star_L < 1/1000 && x_star_R > 1/1000
        Nu_m = integral(@(x) a1*x.^(-c1)+b1,x_star_L,1/80) + integral(@(x) 2.975+a2*x.^(-b2).*exp(-c2*x),1/80,x_star_U);
    end
    Nu_m = Nu_m/(x_star_U-x_star_L);
end