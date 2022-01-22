function d = setup()
    % load data
    S = load('taxi_data\taxi_data.mat');

    d = struct();
    
    % coefficients for objective
    d.p.rho = 0.01;
    d.p.rho_2 = 0.001;
    d.p.rho_c = 0;
    
    % prediction horizon
    d.p.N_NMPC = 20;
    
    % number of vehicles
    d.p.veh = 40; 
    
    % number of nodes
    d.p.nodes = size(S.tij, 1);
    
    % list of nodenames
    d.p.nodenames = [1:d.p.nodes]; 
    
    % length of a single timestep (seconds)
    d.p.T = S.t_step;
    
    % total demand over day (PAX)
    d.p.total_d = sum(S.cij_2, 'all');
    
    % total number of timesteps
    d.p.t_final = size(S.cij_2, 3);

    % rate of discharge for vehicles while moving
    d.p.a_d = double(d.p.T/6)*0.0037;
    
    % rate of charge for vehicles while waiting
    d.p.a_c = 3*d.p.a_d;
    
    % travel time from node i to j
    d.p.tij = S.tij;
    
    % maximum distance between two nodes in timesteps plus 2 for buffer
    d.p.t_max = max(d.p.tij, [], 'all') + 2;
    
    % indicator if travel from node i to j takes T - 1 timesteps
    d.p.tijT = zeros([d.p.nodes, d.p.nodes, d.p.t_max]);
    for i = 1:d.p.nodes
        for j = 1:d.p.nodes
            d.p.tijT(i, j, d.p.tij(i, j) + 1) = 1;
        end
    end
    % shift tijT along T axis by 1 to easily calculate 'tij - 1 = T' later
    d.p.tijT = cat(3, d.p.tijT(:, :, 2:d.p.t_max), zeros([d.p.nodes, d.p.nodes, 1]));
    
    % initial values for dij
    % number of customers demanding travel from node i to node j
    d.p.dij0 = zeros([d.p.nodes, d.p.nodes]);
    
    % initial values for uik
    % indicator if vehicle k stayed in node i at previous step
    d.p.uik0 = zeros([d.p.nodes, d.p.veh]); 
    for i = 1:d.p.nodes
        d.p.uik0(i, floor(d.p.veh/d.p.nodes)*(i - 1) + [1:floor(d.p.veh/d.p.nodes)]) = 1;
    end
    if d.p.nodes*floor(d.p.veh/d.p.nodes) < d.p.veh
        for i = d.p.nodes*floor(d.p.veh/d.p.nodes) + 1:d.p.veh
            d.p.uik0(1, i) = 1;
        end
    end
    
    % initial values for qk
    % initially all set to 1 (full)
    d.p.qk0 = ones([d.p.veh, 1]);
    
    % initial values for pikT
    % indicator if vehicle k is T steps away from node i
    d.p.pikT0 = zeros([d.p.nodes, d.p.veh, d.p.t_max]);
    
    % pre-allocate memory
    d.s.dij = NaN(d.p.nodes, d.p.nodes, d.p.t_final);
    d.s.pikT = NaN(d.p.nodes, d.p.veh, d.p.t_max, d.p.t_final);
    d.s.uik = NaN(d.p.nodes, d.p.veh, d.p.t_final);
    d.s.vijk = NaN(d.p.nodes, d.p.nodes, d.p.veh, d.p.t_final);
    d.s.wijk = NaN(d.p.nodes, d.p.nodes, d.p.veh, d.p.t_final);
    d.s.qk = NaN(d.p.veh, d.p.t_final);
    
    % travel demand
    d.p.cij = zeros(d.p.nodes, d.p.nodes, d.p.t_final + d.p.N_NMPC);
    d.p.cij(:, :, 1:d.p.t_final) = S.cij_2(:, :, 1:d.p.t_final);
    
    % set initial state
    d.s.dij(:, :, 1) = d.p.dij0;
    d.s.pikT(:, :, :, 1) = d.p.pikT0;
    d.s.uik(:, :, 1) = d.p.uik0;
    d.s.qk(:, 1) = d.p.qk0;
   
    % state constraints
    d.p.dij_min = 0;
    d.p.pikT_min = 0;
    d.p.pikT_max = 1;
    d.p.uik_min = 0;
    d.p.uik_max = 1;
    
    % control input constraints
    d.p.vijk_min = 0;
    d.p.vijk_max = 1;
    d.p.wijk_min = 0;
    d.p.wijk_max = 1;
    
end

