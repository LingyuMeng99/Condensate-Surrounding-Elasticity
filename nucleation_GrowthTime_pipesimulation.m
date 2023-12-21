rng('shuffle');

% We only record the condensate growth and dissolution time instead of
% volume with time to decrease the memory load. 
% Some parameters are loaded outside this paragram. 

gamma = 0.1; % Surface tension

c0 = 1; % Saturation concentration, which is 0 by definition
cout = 2; % Initial outside concentration
cin = 10; % Inside concentration

% Initial value
g_ini = cin/c0*log(cout/c0) - (cout/c0-1); % Initial nucleation affinity
Rc = 2 * gamma ./ g_ini;

Rn = kRn * Rc * (ones(N,1)+noise_ratio*(2*rand(N,1)-1)); % Nucleation radius for each condensate
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
tnext = 0;
count = 1;

t_rec_ini = 2e4;
t_growth_list_sum = cell(N,1); % Growth time record
t_dissolution_list_sum = cell(N,1); % Dissolution time record
V_final_sum = cell(N,1); % Final volume record
V_intgrowth_sum = cell(N,1); % Integral of volume with time during growth record
V_int_list = zeros(N,1); % Temporary integral list
check_int_on = ones(N,1) > 0; % 1 for integral on
growth_num_sum = zeros(N,1); % Total times of growth
V_before1 = NaN*ones(N,1); % Volume at t-1
V_before2 = NaN*ones(N,1); % Volume at t-2

dt = 1;
while max(growth_num_sum) <= growth_num_max
    
    if t >= tnext
        if t > t_rec_ini

            % Find dissolved condensates
            dissolution_idx = V_before1 > 1.001*Vn  &  abs(V-Vn) < 0.0001*Vn;
            for idx = find(dissolution_idx ~= 0)'
                t_dissolution_list_sum{idx} = [t_dissolution_list_sum{idx},t];
                V_final_sum{idx} = [V_final_sum{idx},V_before1(idx)];
                V_intgrowth_sum{idx} = [V_intgrowth_sum{idx},V_int_list(idx)];
                check_int_on(idx) = 0; % Integral stop for dissolved condensates
                V_int_list(idx) = 0;
            end
            
            % Find growth initiation
            growth_idx = V_before1 < 1.0001*Vn  &  V-V_before1 > 0.0001*Vn  &  V_before1-V_before2 <= 0.0001*Vn;
            for idx = find(growth_idx~=0)'
                t_growth_list_sum{idx} = [t_growth_list_sum{idx},t-1];
                check_int_on(idx) = 1; % Integral start after growth
            end
            growth_num_sum(growth_idx) = growth_num_sum(growth_idx)+1;

            V_int_list(check_int_on) = V_int_list(check_int_on) + V(check_int_on);
        end

        V_before2 = V_before1;
        V_before1 = V;

        tnext = tnext + dt_rec;
        count = count + 1;

    end

    cout = (N_tot-cin*(sum(V)-Vn_tot))/(V_tot-sum(V));
    g = cin/c0*log(cout/c0) - (cout/c0-1);    
    R = (3*V/(4*pi)).^(1/3);
    Pc = 2*gamma./R + E + Y.*R;
    
    dV = dt * R .* (g-Pc);
    dis_judge = rand(N,1) < k_dis.*V*dt; % Decide whether a condensate is dissolved
    dV(dis_judge) = Vn(dis_judge) - V(dis_judge);
    
    t = t + dt;
    V = V + dV;
    V(V<=Vn) = Vn(V<=Vn);
        
end


