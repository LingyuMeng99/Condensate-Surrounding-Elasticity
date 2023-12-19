rng('shuffle');

% Simulation parameters

N = 10000;
gamma = 0.1; % Surface tension

c0 = 1; % Saturation concentration, which is 0 by definition
cout = 2; % Initial outside concentration
cin = 10; % Inside concentration

% % Power-law random E
E_max = 2;
theta = 0; %x^theta distribution
E = E_max*rand(N,1).^(1/(theta+1));
E = sort(E);

% % Truncated Gaussian random E
% E_mean = 2;
% E_range = E_mean;
% E_var = 1;
% E = zeros(N,1);
% for i = 1:N
%     E_pick = normrnd(E_mean,E_var);
%     while E_pick < E_mean-E_range
%         E_pick = normrnd(E_mean,E_var);
%     end
%     E(i) = E_pick;
% end
% min(E)

% % Exponential random E
% E_mean = 1;
% E_range = E_mean;
% E = (E_mean-E_range) + exprnd(E_range,[N,1]);

% % -3 power-law random E
% E_mean = 1;
% E_range = E_mean;
% E = E_mean - 2*E_range + E_range*(1-rand(N,1)).^(-1/(3-1));

% % Lognormal random E
% sigma_lgn = 0.4;
% mu_lgn = -1.6;
% E = exp(mu_lgn + sigma_lgn.*randn(N,1));

% Beyond neo-Hookean confining pressure
Y = 0.0;

% Condensate dissolution rate per unit volume
k_dis = 0;


% Initial value
g_ini = cin/c0*log(cout/c0) - (cout/c0-1); % Initial nucleation affinity
Rc = 2 * gamma ./ g_ini;

Rn = 3 * Rc * (ones(N,1)+0*(2*rand(N,1)-1)); % Nucleation radius for each condensate
Vn = 4/3*pi * Rn.^3;
Vn_tot = sum(Vn); % Total nucleation volume

V_tot = 1000 * N; % Total volume
N_tot = cout*(V_tot - Vn_tot);  % Total number of molecules


% Simulation 

g = g_ini;
R = Rn;
V = 4/3*pi*R.^3;
Pc = 2*gamma./R + E + Y.*R;

t = 0;
t_total = 1e4;
dt_rec = 1;
total_rec = round(t_total/dt_rec);
V_time = zeros(N,total_rec); % Record volume with time
R_time = zeros(N,total_rec); % Record radius with time
ct = []; % Record average values
tnext = 0;
count = 1;

dt = 1;
while t <= t_total
    
    if t >= tnext
        shrink = ((V-Vn)<=0.001*Vn);
        Rc = 2*gamma/g;

        R_tmean = sum(R.*(1-shrink))/sum(1-shrink); % Direct average R
        R_tmeanV = sum(R.*V)/sum(V); % Volume-weighted average R
        R_tvar = sum(R.^2.*(1-shrink))/sum(1-shrink)-R_tmean^2; % Volume-weighted variation of R
        E_tmean = sum(E.*(1-shrink))/sum(1-shrink); % Direct average E
        E_tmeanV = sum(E.*V)/sum(V); % Volume-weighted average E
        E_tvar = sum(E.^2.*V)/sum(V)-E_tmeanV.^2; % Volume-weighted variation of E
        ER3_tmeanV = sum((E-min(E)).^(1+theta).*R.^3.*V)/sum(V); % Conservation validation

        
        ct = [ct;[t cout g Rc R_tmean R_tmeanV R_tvar E_tmean E_tmeanV E_tvar ER3_tmeanV]];
        V_time(:,count) = V;
        R_time(:,count) = R;

        tnext = tnext + dt_rec;
        count = count + 1;
    end

    cout = (N_tot-cin*(sum(V)-Vn_tot))/(V_tot-sum(V)); % Outside concentration
    g = cin/c0*log(cout/c0) - (cout/c0-1); % Nucleation affinity
    R = (3*V/(4*pi)).^(1/3);
    Pc = 2*gamma./R + E + Y.*R; % Confining pressure
    
    dV = dt * R .* (g-Pc);
    dis_judge = rand(N,1) < k_dis.*V; % Decide whether a condensate is dissolved
    dV(dis_judge) = Vn(dis_judge) - V(dis_judge); % Dissolved condensate V = Vn
    
    t = t + dt;
    V = V + dV;
    V(V<=Vn) = Vn(V<=Vn);

end

%% Plot R-t

loglog(ct(:,1),ct(:,6),'Linewidth',1.5)
hold on
loglog(ct(:,1),(0.96*5*ct(:,1)).^(1/5),'--','Linewidth',2,'Color',[0.6 0.6 0.6])
loglog(ct(:,1),(0.018*ct(:,1)).^(1/3),'--','Linewidth',2,'Color',[0.6 0.6 0.6])
set(gca,'Fontsize',15)
xlabel('$t$','Fontsize',20,'Interpreter','latex')
ylabel('$\langle R \rangle$','Fontsize',20,'Interpreter','latex')
% text(3e4,5,'$t^{1/3}$','Fontsize',20,'Interpreter','latex')


