clear all
close all
clc
eps_steel = 2*10^(-3)/1000;
D = 0.1;
Re_x = 10^4;
f_darcy = 0.001;
L = 10;

N_seg = 200;
L_seg = L/N_seg;
x_edges = L_seg*(0:N_seg);
x_centers = (x_edges(1:end-1)+x_edges(2:end))/2;
Pr = 0.71;
Nu_inf_turb = (f_darcy/8).*(Re_x-1000).*Pr./(1+12.7.*(f_darcy/8).^0.5.*(Pr.^(2/3)-1)); % Eq 7.41 from lienhard/7.76 from Sekulic 

for ii = 1:N_seg
    x_L = x_edges(ii);
    x_R = x_edges(ii+1);
    x_c = x_centers(ii);

    if ii == 1
        Nu_x_turb_1(ii) = Nu_inf_turb*(1 + 0.9757/(x_R/D)^0.760);
        Nu_x_turb_2(ii) = Nu_inf_turb*(1 + (x_R/D)^(-0.9)/Pr^(1/6)*(0.68+3000/Re_x^0.81));
    else
        Nu_x_turb_1(ii) =   1/L_seg*Nu_inf_turb*( x_R*(1 + 0.9757/(x_R/D)^0.760) ...
                                                 -x_L*(1 + 0.9757/(x_L/D)^0.760)); % Average value from x_L to x_R (Lienhard)
        Nu_x_turb_2(ii) = 1/L_seg*Nu_inf_turb*(  x_R*(1 + (x_R/D)^(-0.9)/Pr^(1/6)*(0.68+3000/Re_x^0.81)) ...
                                                -x_L*(1 + (x_L/D)^(-0.9)/Pr^(1/6)*(0.68+3000/Re_x^0.81))); % Average value from x_L to x_R
    end
    Nu_mean_turb_2(ii) = (1 + (x_R/D)^(-0.9)/Pr^(1/6)*(0.68+3000/Re_x^0.81))*Nu_inf_turb;
    %Nu_mean_turb_vpa(ii) = 1/x_R*double(vpaintegral(get_Nu_sym,x,[0,x_R])); % Use mean value for the first cell -> accounts for asymptotic behaviour.
    
end
disp("End")
Nu_x_mean_1 = cumsum(Nu_x_turb_1*L_seg)./x_edges(2:end);
Nu_x_mean_2 = cumsum(Nu_x_turb_2*L_seg)./x_edges(2:end);
f = figure;
hold on
p1 = plot(x_edges(2:end),Nu_x_mean_1,'-','Color',"#0072BD");
p2 = plot(x_edges(2:end),Nu_x_mean_2,'-','Color',"#D95319");
%plot(x_edges(2:end),Nu_mean_turb,':','LineWidth',1.5,'Color',	"#7E2F8E");
p3 = plot(x_centers,Nu_x_turb_1,'--','Color',"#0072BD");
p4 = plot(x_centers,Nu_x_turb_2,'--','Color',"#D95319");
legend([p1,p2,p3,p4],"Mean Lienhard","Mean Sekulic","Local Lienhard","Local Sekulic")
xlim([0.1,gca().XLim(2)])
ylim([0.8,1.5])
