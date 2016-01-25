% for a fixed system size, this script loops over different
% levels (location of the marked vertex) and finds the optimal
% value of the search parameter according to two different
% defintions:
% best_gamma_p: gamma that maximizes probability at the m.v.
% best_gamma_tp: gamma that minimizes t/p


num_levels = 16
max_level = 8;
levels = [1:max_level];
%levels = [ 1 num_levels/4 num_levels/2 3*num_levels/4 ];
%max_level = length(levels);

gammas = [[0.1:0.1:0.3],[0.4:0.05:1.2],[1.4:0.2:4.4]];
%gammas = [0.5:0.01:1.1];
%gammas = [[0.65:0.001:0.68],[1]];

num_times = 100;	
min_time = 3; % This value was chosen based on numerical experiments
max_time = sqrt(2^num_levels);

best_gamma_p = zeros(max_level,1);
best_gamma_tp = zeros(max_level,1);


for l = levels

	tic;
    
    	max_time = 5*max_time+8;
	time_steps =min(max(2,(max_time-min_time)/num_times),100);    	
	times = [ min_time:time_steps:max_time ];
   	%times = max_time.*rand(1,num_times);

	[ Lb, imv, marked_state, initial_state ] = generate_reduced_L_directly( num_levels, l );
	Lb = sparse(Lb);
	marked_state = sparse(marked_state);

	max_p_values = zeros(1,length(gammas));
	min_tp_values = zeros(1,length(gammas));

	gam_count = 1;
	maxtime_idx = 0;	

	for gam = gammas

		Hb=gam.*Lb-diag(marked_state);
		
		p_values = zeros(1,length(times));	
		tp_values = zeros(1,length(times));	
		
		t_count = 1;

		for t = times
			
			prob = abs( dot( expm(-i*Hb*t) * initial_state, marked_state ) )^2;
			p_values(1,t_count) = prob;
            		tp_values(1,t_count) = t/prob;
			
		        t_count = t_count + 1;

		end

		max_p_values( gam_count ) = max( p_values );
		
		[mn,idx] = min( tp_values );
		min_tp_values( gam_count ) = mn;

		%fprintf('\t%d %d\n',length(times),idx)

		maxtime_idx = max(maxtime_idx,idx-1);
		
		gam_count = gam_count + 1;

	end

	[mx,idx] = max( max_p_values);
	best_gamma_p(l) = gammas(idx);

	[mx,idx] = min( min_tp_values);
	best_gamma_tp(l) = gammas(idx);

    	max_time=times(maxtime_idx+1);        		

    	t0=toc;
	fprintf('level: %d  elapsed time: %.2f  maxtime_idx: %d  max_time: %.2f \n', l, t0, maxtime_idx, max_time);
	if maxtime_idx <= 1
		fprintf('\t warning: used trivial times for all gammas \n');
	end        	

  	max_time = max_time/2;

end


plot(levels,best_gamma_p,'r',levels,best_gamma_tp,'b')
xlabel('level l of the marked vertex')
ylabel('gamma')
title(sprintf('gamma_{opt} as function of l for n=%d (...)',num_levels))
