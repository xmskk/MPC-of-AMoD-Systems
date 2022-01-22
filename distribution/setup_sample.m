function d = setup_sample()

    d = struct();
    
    % coefficients for objective
    d.p.rho = 0.01;
    d.p.rho_2 = 0.001;
    d.p.rho_c = 0;
    
    % prediction horizon
    d.p.N_NMPC = 10;
    
    % number of vehicles
    d.p.veh = 30; 
    
    % number of nodes
    d.p.nodes = 10;
    
    % list of nodenames
    d.p.nodenames = [1:d.p.nodes]; 
    
    % length of a single timestep (seconds)
    d.p.T = 6;
    
    % total number of timesteps
    d.p.t_final = 30;
    
    % rate of charge for vehicles while waiting
    d.p.a_c = 0.2;
    
    % rate of discharge for vehicles while moving
    d.p.a_d = 0.1;
    
    % travel time from node i to j
    d.p.tij = [0 2 1 3 4 4 5 6 5 7;...
        2 0 3 1 2 3 3 4 6 5;...
        1 3 0 2 3 2 4 5 4 6;...
        3 1 2 0 1 2 2 3 4 4;...
        4 2 3 1 0 2 1 2 5 3;...
        4 3 2 2 2 0 2 3 3 4;...
        5 3 4 2 1 2 0 1 4 2;...
        6 4 5 3 2 3 1 0 5 2;...
        5 6 4 4 5 3 4 5 0 5;...
        7 5 6 4 3 4 2 2 5 0]; 
    
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
    for i = d.p.nodenames
        n = randsample([15:d.p.veh], 1);
        fprintf('%d users at node %d\n', n, i)
        for m = 1:n
            j = randsample(d.p.nodenames(d.p.nodenames ~= i), 1);
            d.p.dij0(i, j) = d.p.dij0(i, j) +  1;
        end
    end
    
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

