function d = nmpc()

    addpath(genpath('yalmip location'))
    addpath(genpath('gurobi location'))
    addpath(genpath(pwd))
    
    d = setup();
    disp('Finished initial setup!')
    
    d = set_con(d, 1);
    disp('Finished setting conditions!')
    
    for t = 1:d.p.t_final
        fprintf('Currently at step #%.0f', t)
        
        d = solve_NMPC(d, t);
        
    end
    
    d = rmfield(d, 'c');
    
    save('results.mat', 'd');
    
end

function d = set_con(d, t)

    % define states as decision variables
    d.c.dij = intvar(d.p.nodes, d.p.nodes, d.p.N_NMPC, 'full');
    d.c.pikT = binvar(d.p.nodes, d.p.veh, d.p.t_max, d.p.N_NMPC, 'full');
    d.c.uik = binvar(d.p.nodes, d.p.veh, d.p.N_NMPC, 'full');
    d.c.qk = sdpvar(d.p.veh, d.p.N_NMPC, 'full');
    
    % define control inputs as decision variables
    d.c.vijk = binvar(d.p.nodes, d.p.nodes, d.p.veh, d.p.N_NMPC - 1, 'full');
    d.c.wijk = binvar(d.p.nodes, d.p.nodes, d.p.veh, d.p.N_NMPC - 1, 'full');
    
    % initialize objective function
    d.c.obj = 0;
    
    % initialize constraints
    d.c.con = [];
    
    % initial state constraint
    d.c.con = [d.c.con, d.c.dij(:, :, 1) == d.s.dij(:, :, t)]; % 1
    d.c.con = [d.c.con, d.c.pikT(:, :, :, 1) == d.s.pikT(:, :, :, t)]; % 2
    d.c.con = [d.c.con, d.c.uik(:, :, 1) == d.s.uik(:, :, t)]; % 3
    d.c.con = [d.c.con, d.c.qk(:, 1) == d.s.qk(:, t)]; % 4
    
    % state variable constraint
    d.c.con = [d.c.con, d.c.dij(:, :, 2:end) >= d.p.dij_min*ones(d.p.nodes, d.p.nodes, d.p.N_NMPC - 1)]; % 5
    d.c.con = [d.c.con, zeros(d.p.veh, d.p.N_NMPC - 1) <= d.c.qk(:, 2:end)]; % 6
    d.c.con = [d.c.con, d.c.qk(:, 2:end) <= ones(d.p.veh, d.p.N_NMPC - 1)]; % 7
        
    % control constraint
    d.c.con = [d.c.con, squeeze(sum(d.c.pikT(:, :, :, 2:end), [1, 3])) <= ones(d.p.veh, d.p.N_NMPC - 1)]; % 8
    d.c.con = [d.c.con, squeeze(sum(d.c.uik(:, :, 2: end), 1)) + squeeze(sum(d.c.pikT(:, :, :, 2:end), [1, 3])) == 1]; % 9
    d.c.con = [d.c.con, squeeze(sum(d.c.uik(:, :, 2: end), 1)) + squeeze(sum(d.c.vijk + d.c.wijk, [1, 2])) <= ones(d.p.veh, d.p.N_NMPC - 1)]; % 10
    d.c.con = [d.c.con, squeeze(sum(d.c.vijk, 3)) <= d.c.dij(:, :, 1:end - 1) + d.p.cij(:, :, t:t + d.p.N_NMPC - 2)]; % 11
    d.c.con = [d.c.con, d.c.qk(:, 2:end) >= d.p.a_d*squeeze(sum(d.c.vijk.*repmat(d.p.tij, 1, 1, d.p.veh, d.p.N_NMPC - 1), [1, 2]))]; % 12
    d.c.con = [d.c.con, d.c.qk(:, 2:end) >= d.p.a_d*squeeze(sum(d.c.wijk.*repmat(d.p.tij, 1, 1, d.p.veh, d.p.N_NMPC - 1), [1, 2]))]; % 13
    
    % evolve dynamic contraint
    d.c.con = [d.c.con, d.c.dij(:, :, 2:end) == d.c.dij(:, :, 1:end - 1) + d.p.cij(:, :, t:t + d.p.N_NMPC - 2) - squeeze(sum(d.c.vijk, 3))]; % 14
    d.c.con = [d.c.con, d.c.pikT(:, :, :, 2:end) == cat(3, d.c.pikT(:, :, 2:d.p.t_max, 1:end - 1), repmat(zeros([d.p.nodes, d.p.veh, 1]), 1, 1, 1, d.p.N_NMPC - 1))...
        + squeeze(sum(permute(repmat(d.c.vijk, 1, 1, 1, 1, d.p.t_max), [1, 2, 3, 5, 4]).*permute(repmat(d.p.tijT, 1, 1, 1, d.p.N_NMPC - 1, d.p.veh), [1, 2, 5, 3, 4]), 1))...
        + squeeze(sum(permute(repmat(d.c.wijk, 1, 1, 1, 1, d.p.t_max), [1, 2, 3, 5, 4]).*permute(repmat(d.p.tijT, 1, 1, 1, d.p.N_NMPC - 1, d.p.veh), [1, 2, 5, 3, 4]), 1))]; % 15
    d.c.con = [d.c.con, d.c.uik(:, :, 2:end) == d.c.uik(:, :, 1:end - 1) + squeeze(d.c.pikT(:, :, 1, 1:end - 1)) - squeeze(sum(d.c.vijk + d.c.wijk, 2))]; % 16
    d.c.con = [d.c.con, d.c.qk(:, 2:end) <= d.c.qk(:, 1:end - 1) + squeeze(d.p.a_c*sum(d.c.uik(:, :, 2:end), 1)) - squeeze(d.p.a_d*sum(d.c.pikT(:, :, :, 2:end), [1, 3]))]; % 17
    d.c.con = [d.c.con, d.c.qk(:, 2:end) <= ones([d.p.veh, d.p.N_NMPC - 1]) - squeeze(d.p.a_d*sum(d.c.pikT(:, :, :, 2:end), [1, 3]))]; % 18
    
    % stage cost
    d.c.obj = d.c.obj + sum(d.c.dij(:, :, 2:end), 'all') + d.p.rho*sum(repmat(d.p.tij, 1, 1, d.p.veh, d.p.N_NMPC - 1).*d.c.wijk, 'all')...
        - d.p.rho_2*sum(d.c.qk(:, 2:end), 'all');
    
    % terminal cost
    d.c.obj = d.c.obj - d.p.rho_c*sum(d.c.qk(:, end), 'all');
    
    % define optimization settings
    d.c.ops = sdpsettings;
    d.c.ops.verbose = 0; % configure level of info solver displays
    d.c.ops.solver = 'gurobi'; % choose solver
    d.c.ops.debug = 1;
    
end

function d = solve_NMPC(d, t)
    tic_c = tic;
    
    % renew constraints that varies per iteration
    % initial state constraint
    d.c.con(1) = (d.c.dij(:, :, 1) == d.s.dij(:, :, t)); % 1
    d.c.con(2) = (d.c.pikT(:, :, :, 1) == d.s.pikT(:, :, :, t)); % 2
    d.c.con(3) = (d.c.uik(:, :, 1) == d.s.uik(:, :, t)); % 3
    d.c.con(4) = (d.c.qk(:, 1) == d.s.qk(:, t)); % 4
    
    % control constraint
    d.c.con(11) = (squeeze(sum(d.c.vijk, 3)) <= d.c.dij(:, :, 1:end - 1) + d.p.cij(:, :, t:t + d.p.N_NMPC - 2)); % 11
    
    % evolve dynamic contraint
    d.c.con(14) = (d.c.dij(:, :, 2:end) == d.c.dij(:, :, 1:end - 1) + d.p.cij(:, :, t:t + d.p.N_NMPC - 2) - squeeze(sum(d.c.vijk, 3))); % 14
    
    % solve NMPC problem and record diagnostics
    d.s.diagnostics{t} = optimize(d.c.con, d.c.obj, d.c.ops);
    
    % record primal and dual constraint violations
    [primalfeas, dualfeas] = check(d.c.con);
    d.s.primalfeas{t} = primalfeas;
    d.s.dualfeas{t} = dualfeas;
    
    % record CPU time
    d.s.CPU_time(t, 1) = toc(tic_c);
    disp(d.s.CPU_time(t, 1))
    
    % extract the control input trajectory
    % from the solution of MPC optimization problem
    vijk_t = value(d.c.vijk);
    wijk_t = value(d.c.wijk);
    
    % extract state variable trajectory
    % from the solution of MPC optimization problem
    dij_t = value(d.c.dij);
    pikT_t = value(d.c.pikT);
    uik_t = value(d.c.uik);
    qk_t = value(d.c.qk);
    
    % assign first element of the solution to the NMPC
    % problem as the control input at time t
    d.s.vijk(:, :, :, t) = vijk_t(:, :, :, 1);
    d.s.wijk(:, :, :, t) = wijk_t(:, :, :, 1);
    
    % assign second element of the solution to the NMPC
    % problem as the state variable at time t
    d.s.dij(:, :, t + 1) = dij_t(:, :, 2);
    d.s.pikT(:, :, :, t + 1) = pikT_t(:, :, :, 2);
    d.s.uik(:, :, t + 1) = uik_t(:, :, 2);
    d.s.qk(:, t + 1) = qk_t(:, 2);
end