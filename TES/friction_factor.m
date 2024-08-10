close all
clear
clc
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaulttextInterpreter','latex')


N_jj = 100000;
all_rel_roughness = [0,0.001,0.002,0.005,0.01,0.02,0.03,0.04,0.05];
N_ii = length(all_rel_roughness);
all_Re = 10.^linspace(3,6,N_jj);

f = figure;
hold on
for ii = 1:N_ii
    rel_roughness = all_rel_roughness(ii);
    for jj = 1:N_jj
        Re = all_Re(jj);
        f_darcy(jj) = GetDarcy(Re,rel_roughness);
        f_darcy_Haaland(jj) = (-1.8*log10((rel_roughness/3.7).^1.11+6.9./Re)).^(-2); % Haaland correlation
        if Re<4000
            f_darcy_lam(jj) = 64/Re;
        end
    end
    plot(all_Re,f_darcy_Haaland,':','Color','k')
    plot(all_Re(all_Re<4000),f_darcy_lam,':','Color','k')
    plot(all_Re,f_darcy,'Color','k')
    text(all_Re(end)^1.01,f_darcy(end),sprintf('%0.3f',rel_roughness));
    all_f_darcy(ii) = f_darcy(end);
end
f_darcy_avg = (all_f_darcy(1)+all_f_darcy(end))/2;
t1 = text(all_Re(end)^1.009,f_darcy_avg,'Relative roughness $\frac{e}{D}$', ...
    'HorizontalAlignment','center','VerticalAlignment','middle');
set(t1,'Rotation',90);
t1.Units = 'pixels';
set(gca,'XScale','log')
xline(2100)
ylabel('Darcy friction factor $f_D$')
xlabel('Reynolds number')
ax = gca;
ax.Units = 'pixels';
f.Position(3) = ax.Position(1)+t1.Position(1)+10;
t1.Units = 'data';

text(10^3.05,0.018,'$f_D = 64/\mathrm{Re}$')
