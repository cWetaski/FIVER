clear
close all
clc

T_air = 0; % K
Ro = 2;
Bi_i = 1;
N_check = 15;
Bi_o = 0;
Ro_all = linspace(1,10,N_check+1);%10.^(linspace(-3,3,N_check));
Ro_all(1) = [];
N_vals = 10000;
x_vals = linspace(0,100,N_vals+1);
x_vals = x_vals(2:end); % Get rid of 0 value
dx = x_vals(2)-x_vals(1); % 'sampling period'
Fs = 1/dx;
peak_freqs = cell(1);
N_peaks = 0;
for nn = 1:N_check
    Ro = Ro_all(nn);
    y_vals = lambda_fun(x_vals,Bi_i,Bi_o,Ro);
    Y_vals = fft(y_vals); % Fast fourier transform of y
    Y_vals_shift = abs(fftshift(Y_vals));
    freqs = Fs/N_vals*(-N_vals/2:N_vals/2-1);
    Y_vals_shift = Y_vals_shift(freqs>=0);
    freqs = freqs(freqs>=0);
    peak_inds = islocalmax(Y_vals_shift);
    peak_vals = Y_vals_shift(peak_inds);
    peak_inds = find(peak_inds & (Y_vals_shift>0.2*max(peak_vals)));
    peak_vals = Y_vals_shift(peak_inds);
    peak_freqs_cur = freqs(peak_inds); 
    N_peaks_cur = length(peak_freqs_cur);
    if true && nn == round(N_check/2)
        t = tiledlayout(1,2);
        nexttile()
        plot(x_vals(x_vals<(10./max(peak_freqs_cur))),y_vals(x_vals<(10./max(peak_freqs_cur))));
        nexttile()
        plot(freqs(freqs<5*max(peak_freqs_cur)),Y_vals_shift(freqs<5*max(peak_freqs_cur)))
    end

    if N_peaks_cur > N_peaks
        N_peaks = N_peaks_cur;
        peak_freqs{N_peaks_cur} = zeros(N_check,1);
        x_vals_saved = x_vals;
        y_vals_saved = y_vals;
        Bi_i_saved = Bi_i;
        Bi_o_saved = Bi_o;
        Ro_saved = Ro;
    end
    for ii = 1:N_peaks_cur
        peak_freqs{ii}(nn) = peak_freqs_cur(ii);
    end
end

f = figure;
hold on
for ii = 1:N_peaks
    plot(Ro_all,peak_freqs{ii})
    p = polyfit(Ro_all,peak_freqs{ii},1);
    if ii == 1
    text(1+(max(Ro_all)-1)*0.6,max(peak_freqs{ii}*0.2),sprintf('y = %0.3fx + %0.3f',p(1),p(2)))
    end
end

%% PART 2 - Checking Solutions
f2 = figure;
N_Fo = 50;
N_rr = 25;
Fo_vals = linspace(0,5,N_Fo+1);
Ro = N_rr*2+1;
Bi_i = 1;
Bi_o = 0;
k = 0.001;
ri = 1;
ro = ri*Ro;
hi = Bi_i*k/ri;
ho = Bi_o*k/ri;
N_lambda = 10;
T_1 = 1;
T_air = 0;

%% Analytic solution
dr = (Ro-1)/N_rr;
R_vals_analytic = [1,(1+dr/2):dr:(Ro-dr/2),Ro];
r_vals_analytic = R_vals_analytic*ri;
lambdas = get_lambdas(N_lambda,Bi_i,Bi_o,Ro);
lambdas2 = get_lambdas2(N_lambda,k,hi,ho,ri,ro);
lambdas3 = get_lambdas3(N_lambda,k,hi,ho,ri,ro);
lambdas4 = get_lambdas4(N_lambda,k,-hi,ho,ri,ro);
x_vals = linspace(0,max(lambdas/ri)*1.25,N_vals+1);
x_vals(1) = [];

y_vals = lambda_fun(x_vals*ri,Bi_i,Bi_o,Ro);
y_vals2 = lambda_fun2(x_vals,k,hi,ho,ri,ro);
y_vals3 = lambda_fun3(x_vals,k,hi,ho,ri,ro);
y_vals4 = lambda_fun4(x_vals,k,-hi,ho,ri,ro);
xline(0,':')
hold on
plot(x_vals,y_vals)
xline(lambdas/ri)
plot(x_vals,y_vals2,'--');
xline(lambdas2,'--')
plot(x_vals,y_vals3,':',LineWidth=2);
xline(lambdas3,':')
plot(x_vals,y_vals4,'-.');
xline(lambdas4,'-.')


W = get_W(R_vals_analytic,lambdas,Bi_i);
W2 = get_W2(r_vals_analytic,lambdas2,k,hi,ri);
E = get_E(lambdas,Bi_i,Bi_o,T_1,T_air,Ro);
E2 = get_E2(lambdas2,k,hi,ho,T_1,T_air,ri,ro);

v_air = T_air-T_1;
ci = -hi*v_air;
co = ho*v_air;
C0 = get_C0(r_vals_analytic,lambdas3,k,hi,ri);
F = get_F(lambdas3,k,hi,ho,ri,ro);
G = get_G(lambdas3,k,hi,ho,ci,co,ri,ro);
H = get_H(r_vals_analytic,k,hi,ho,ci,co,ri,ro);

C04 = get_C04(r_vals_analytic,lambdas4,k,-hi,ri);
F4 = get_F4(lambdas4,k,-hi,ho,ri,ro);
G4 = get_G4(lambdas4,k,-hi,ho,ci,co,ri,ro);


f = figure;
for ii = 1:N_Fo
    sol = HollowCylinderConvectionTransient(Fo_vals(ii)*Ro^2,lambdas,W,E);
    sol2 = HollowCylinderConvectionTransient2(Fo_vals(ii)*ro^2,lambdas2,W2,E2);
    v3 = HollowCylinderConvectionTransient3(Fo_vals(ii)*ro^2,lambdas3,C0,F,G,H);
    v4 = HollowCylinderConvectionTransient4(Fo_vals(ii)*ro^2,lambdas4,C04,F4,G4,H);
    T_R_analytic1 = T_air + sol;
    T_R_analytic2 = T_air - sol2;
    T_R_analytic3 = T_1 + v3;
    T_R_analytic4 = T_1 + v4;
    plot(R_vals_analytic,T_R_analytic1);
    hold on
    plot(R_vals_analytic,T_R_analytic2,'--');
    plot(R_vals_analytic,T_R_analytic3,':');
    plot(R_vals_analytic,T_R_analytic4,'-.');
    text(1+max(R_vals_analytic)*0.2,max(T_1,T_air)*0.1,sprintf('Fo = %0.2f',Fo_vals(ii)))
    ylim([min(T_1,T_air),max(T_1,T_air)])
    hold off
end

function sol = HollowCylinderConvectionTransient(Fo,lambdas,W,E)
        sol = sum(E.*exp(-(lambdas.^2).*Fo).*W,1);
end

function sol = HollowCylinderConvectionTransient2(Fo,lambdas,W,E)
        sol = sum(E.*exp(-(lambdas.^2).*Fo).*W,1);
end

function v = HollowCylinderConvectionTransient3(Fo,lambdas,C0,F,G,H)
    v = H-pi*sum(exp(-lambdas.^2*Fo).*C0.*F.^(-1).*G,1);
end

function v = HollowCylinderConvectionTransient4(Fo,lambdas,C0,F,G,H)
    v = H-pi*sum(exp(-lambdas.^2*Fo).*F.*C0.*G,1); 
end

function val = lambda_fun(lambda,Bi_i,Bi_o,Ro)
 val = (Bi_i*J0(lambda)+lambda.*J1(lambda)).*(Bi_o*Y0(lambda*Ro)-lambda.*Y1(lambda*Ro)) ...
     - (Bi_o*J0(lambda*Ro)-lambda.*J1(lambda*Ro)).*(Bi_i*Y0(lambda)+lambda.*Y1(lambda));
end

function val = lambda_fun2(lambda,k,hi,ho,ri,ro)
    val = ((hi/k)*J0(lambda*ri)+lambda.*J1(lambda*ri)).*((ho/k)*Y0(lambda*ro)-lambda.*Y1(lambda*ro))...
         -((ho/k)*J0(lambda*ro)-lambda.*J1(lambda*ro)).*((hi/k)*Y0(lambda*ri)+lambda.*Y1(lambda*ri));
end

function lambdas = get_lambdas(N,Bi_i,Bi_o,Ro)
    syms lambda
    eqn = (Bi_i*J0(lambda)+lambda*J1(lambda))*(Bi_o*Y0(lambda*Ro)-lambda*Y1(lambda*Ro))...
        - (Bi_o*J0(lambda*Ro)-lambda*J1(lambda*Ro))*(Bi_i*Y0(lambda)+lambda*Y1(lambda)) ==0;
    safety_factor = 10;
    sampling_freq = 0.1592*(Ro-1)*safety_factor; % Determined this value of 0.1592 seperately
    sampling_period = 1/sampling_freq;
    lambdas = zeros(N,1);
    n = 0;
    lambda_found = 0;
    t_prev = 1e-6; % offset it a bit
    val_prev = lambda_fun(t_prev,Bi_i,Bi_o,Ro); % First value very close to 0
    while lambda_found < N
        n = n+1;
        t_cur = t_prev+sampling_period;
        val_cur = lambda_fun(t_cur,Bi_i,Bi_o,Ro);
        if val_cur*val_prev<0
            lambda_cur = double(vpasolve(eqn,lambda,[t_prev t_cur]));
            err = lambda_fun(lambda_cur,Bi_i,Bi_o,Ro);
            if err<1e-2
                lambda_found = lambda_found +1;
                lambdas(lambda_found) = lambda_cur;
            end
        end
        val_prev = val_cur;
        t_prev = t_cur;
    end
end

function lambdas = get_lambdas2(N,k,hi,ho,ri,ro)
    syms lambda
    eqn = ((hi/k)*J0(lambda*ri)+lambda.*J1(lambda*ri)).*((ho/k)*Y0(lambda*ro)-lambda.*Y1(lambda*ro))...
         -((ho/k)*J0(lambda*ro)-lambda.*J1(lambda*ro)).*((hi/k)*Y0(lambda*ri)+lambda.*Y1(lambda*ri))==0;
    safety_factor = 10; % Used to make the root finding quicker since we search a much smaller region for the root
    sampling_freq = 0.1592*2*(ro-ri)*safety_factor; % Determined this value of 0.16 seperately
    sampling_period = 1/sampling_freq;
    lambdas = zeros(N,1);
    n = 0;
    lambda_found = 0;
    t_prev = 1e-6;
    val_prev = lambda_fun2(t_prev,k,hi,ho,ri,ro); % First value very close to 0
    while lambda_found < N
        n = n+1;
        t_cur = t_prev+sampling_period;
        val_cur = lambda_fun2(t_cur,k,hi,ho,ri,ro);
        if val_cur*val_prev<0
            lambda_cur = double(vpasolve(eqn,lambda,[t_prev t_cur]));
            err = lambda_fun2(lambda_cur,k,hi,ho,ri,ro);
            if err<1e-2
                lambda_found = lambda_found +1;
                lambdas(lambda_found) = lambda_cur;
            end
        end
        val_prev = val_cur;
        t_prev = t_cur;
    end
end

function W = get_W(R,lambdas,Bi_i)
    W = -(Bi_i.*Y0(lambdas) + lambdas.*Y1(lambdas)).*J0(lambdas.*R) ...
        +(lambdas.*J1(lambdas)+Bi_i*J0(lambdas)).*Y0(lambdas.*R);
end

function W = get_W2(r,lambdas,k,hi,ri)
    W = -((hi/k)*Y0(lambdas*ri)+lambdas.*Y1(lambdas*ri)).*J0(lambdas.*r)...
        +(lambdas.*J1(lambdas.*ri)+(hi/k)*J0(lambdas.*ri)).*Y0(lambdas.*r);
end

function E = get_E(lambdas,Bi_i,Bi_o,T_1,T_f,Ro)
    intfun = @(Rint,lambda) Rint.*(T_1-T_f).*get_W(Rint,lambda,Bi_i);
    intvals = zeros(length(lambdas),1);
    for ii = 1:length(lambdas)
        intfun_cur = @(Rint) intfun(Rint,lambdas(ii));
        intvals(ii) = integral(intfun_cur,1,Ro);
    end
    E = pi^2/2*lambdas.^2.*(Bi_o*J0(lambdas*Ro)-lambdas.*J1(lambdas*Ro)).^2.*intvals...
        .*((lambdas.^2+Bi_o^2).*(Bi_i*J0(lambdas)+lambdas.*J1(lambdas)).^2 ...
        -(lambdas.^2+Bi_i^2).*(Bi_o*J0(lambdas*Ro)-lambdas.*J1(lambdas*Ro)).^2).^(-1);    
end

function E = get_E2(lambdas,k,hi,ho,T_1,T_air,ri,ro)
    intfun = @(r,lambda) r.*(T_air-T_1).*get_W2(r,lambda,k,hi,ri);
    intvals = zeros(length(lambdas),1);
    for ii = 1:length(lambdas)
        intfun_cur = @(r) intfun(r,lambdas(ii));
        intvals(ii) = integral(intfun_cur,ri,ro);
    end
    E = pi^2*lambdas.^2/2 ...
    .*((ho/k)*J0(lambdas*ro)-lambdas.*J1(lambdas*ro)).^2.*intvals ...
    .*((lambdas.^2+(ho/k)^2).*(hi/k*J0(lambdas*ri)+lambdas.*J1(lambdas*ri)).^2 ...
      -(lambdas.^2+(hi/k)^2).*(ho/k*J0(lambdas*ro)-lambdas.*J1(lambdas*ro)).^2).^(-1);
end

function val = lambda_fun3(lambda,k,hi,ho,ri,ro)
    val = (k*lambda.*J1(ri*lambda)+hi*J0(ri*lambda)).*(k*lambda.*Y1(ro*lambda)-ho*Y0(ro*lambda)) ...
         -(k*lambda.*Y1(ri*lambda)+hi*Y0(ri*lambda)).*(k*lambda.*J1(ro*lambda)-ho*J0(ro*lambda));
end

function val = lambda_fun4(lambda,k2,k3i,k3o,a,b)
    val = (-k2*lambda.*J1(a*lambda)+k3i*J0(a*lambda)).*(-k2*lambda.*Y1(b*lambda)+k3o*Y0(b*lambda)) ...
         -(-k2*lambda.*Y1(a*lambda)+k3i*Y0(a*lambda)).*(-k2*lambda.*J1(b*lambda)+k3o*J0(b*lambda));
end

function lambdas = get_lambdas3(N,k,hi,ho,ri,ro)
    syms lambda
    eqn = (k*lambda*J1(ri*lambda)+hi*J0(ri*lambda))*(k*lambda*Y1(ro*lambda)-ho*Y0(ro*lambda)) ...
         -(k*lambda*Y1(ri*lambda)+hi*Y0(ri*lambda))*(k*lambda*J1(ro*lambda)-ho*J0(ro*lambda)) == 0;
    safety_factor = 10;
    sampling_freq = 0.16*2*(ro-ri); % Determined this relationship seperately
    sampling_period = 1/sampling_freq;
    sampling_period = sampling_period/safety_factor;
    lambdas = zeros(N,1);
    n = 0;
    lambda_found = 0;
    t_prev = 1e-6;
    val_prev = lambda_fun3(t_prev,k,hi,ho,ri,ro); % First value very close to 0
    while lambda_found < N
        n = n+1;
        t_cur = t_prev+sampling_period;
        val_cur = lambda_fun3(t_cur,k,hi,ho,ri,ro);
        if val_cur*val_prev<0
            lambda_cur = double(vpasolve(eqn,lambda,[t_prev t_cur]));
            err = lambda_fun3(lambda_cur,k,hi,ho,ri,ro);
            if err<1e-2
                lambda_found = lambda_found +1;
                lambdas(lambda_found) = lambda_cur;
            end
        end
        val_prev = val_cur;
        t_prev = t_cur;
    end
end

function lambdas = get_lambdas4(N,k2,k3i,k3o,a,b)
    syms lambda
    eqn = (-k2*lambda*J1(a*lambda)+k3i*J0(a*lambda))*(-k2*lambda*Y1(b*lambda)+k3o*Y0(b*lambda)) ...
         -(-k2*lambda*Y1(a*lambda)+k3i*Y0(a*lambda))*(-k2*lambda*J1(b*lambda)+k3o*J0(b*lambda)) == 0;
    safety_factor = 10;
    sampling_freq = 2*0.16*(b-a); % Determined this value of 0.16 seperately
    sampling_period = 1/sampling_freq;
    sampling_period = sampling_period/safety_factor;
    lambdas = zeros(N,1);
    n = 0;
    lambda_found = 0;
    t_prev = 1e-6;
    val_prev = lambda_fun4(t_prev,k2,k3i,k3o,a,b); % First value very close to 0
    while lambda_found < N
        n = n+1;
        t_cur = t_prev+sampling_period;
        val_cur = lambda_fun4(t_cur,k2,k3i,k3o,a,b);
        if val_cur*val_prev<0
            lambda_cur = double(vpasolve(eqn,lambda,[t_prev t_cur]));
            err = lambda_fun4(lambda_cur,k2,k3i,k3o,a,b);
            if err<1e-2
                lambda_found = lambda_found +1;
                lambdas(lambda_found) = lambda_cur;
            end
        end
        val_prev = val_cur;
        t_prev = t_cur;
    end
end

function C0 = get_C0(r,lambdas,k,hi,ri)
    C0 = J0(r.*lambdas).*(k*lambdas.*Y1(ri*lambdas)+hi*Y0(ri*lambdas))...
        -Y0(r.*lambdas).*(k*lambdas.*J1(ri*lambdas)+hi*J0(ri*lambdas));
end

function C0 = get_C04(r,lambdas,k,k3i,a)
    C0 = J0(r.*lambdas).*(-k*lambdas.*Y1(a*lambdas)+k3i*Y0(a*lambdas))...
        -Y0(r.*lambdas).*(-k*lambdas.*J1(a*lambdas)+k3i*J0(a*lambdas));
end

function F = get_F(lambdas,k,hi,ho,ri,ro)
    denom = (k*lambdas.*J1(ro*lambdas)-ho*J0(ro*lambdas));
    F = (k^2*lambdas.^2+ho^2).*(k*lambdas.*J1(ri*lambdas)+hi*J0(ri*lambdas)).^2 ...
       -(k^2*lambdas.^2+hi^2).*(k*lambdas.*J1(ro*lambdas)-ho*J0(ro*lambdas)).^2;
    F = F./denom;
end

function F = get_F4(lambda,k,k3i,k3o,a,b)
% Note: in this formulation, k3i is -hi, not +hi (the negative of the heat transfer coefficient on the inner surface)
    num = k3o*J0(b*lambda)-k*lambda.*J1(b*lambda);
    denom = (k3o*J0(b*lambda)-k*lambda.*J1(b*lambda)).^2.*(k3i^2+k.^2*lambda.^2) ...
           -(k3i*J0(a*lambda)-k*lambda.*J1(a*lambda)).^2.*(k3o^2+k.^2*lambda.^2);
    F = num./denom;
end


function G = get_G(lambdas,k,hi,ho,ci,co,ri,ro)
    G = ci*(k*lambdas.*J1(ro*lambdas)-ho*J0(ro*lambdas)) ...
       -co*(k*lambdas.*J1(ri*lambdas)+hi*J0(ri*lambdas));
end

function G = get_G4(lambdas,k,k3i,k3o,k4i,k4o,a,b)
% Note: in this formulation, k3i is -hi, not +hi (the negative of the heat transfer coefficient on the inner surface)
    G = k4i*(k3o*J0(b*lambdas)-k*lambdas.*J1(b*lambdas)) ...
       -k4o*(k3i*J0(a*lambdas)-k*lambdas.*J1(a*lambdas));
end

function H = get_H(r,k,hi,ho,ci,co,ri,ro)
    H = (-ri*ci*(k-ro*ho*log(r/ro))+ro*co*(k+ri*hi*log(r/ri))) ...
        /(ri*hi*k+ro*k*ho+ri*ro*hi*ho*log(ro/ri));
end

function J0 = J0(x)
    J0 = besselj(0,x);
end

function J1 = J1(x)
    J1 = besselj(1,x);
end

function Y0 = Y0(x)
    Y0 = bessely(0,x);
end

function Y1 = Y1(x)
    Y1 = bessely(1,x);
end