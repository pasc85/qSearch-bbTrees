%% this script produces a plot of the largest probability concentrations
%% that can by achieved by different values of the search parameter

clear;

n_list = [ 8 12 16 ]; 	% list of system sizes (change plot command at the
			% bottom when changing the length of this list)

gammas = [[0.0:0.05:0.45],[0.5:0.01:0.9],[0.95:0.05:2.0]];	% range of gammas
%gammas = [0.64:0.001:0.69];

num_times = 200;	% number of time samples

% max_p_values saves the data
max_p_values = zeros(length(n_list),length(gammas));
n_count = 1;

for n = n_list	
	
	tic;
	
	% use equally spaced time points or random ones
	max_time = 2^n;		
	times = [ 1:max_time/num_times:max_time ]; 
	%times = max_time.*rand(1,num_times);

	gam_count = 1;	

	l = n/2;	% level of the marked vertex,
	%(change to 1,2,3,n/4,n/2,n-2,n-1,n etc; has to be an integer)

	% generate reduced Hamiltonian and vectors
	[ Lb, imv, marked_state, initial_state ] = generate_reduced_L_directly( n, l );
	Lb = sparse(Lb);
	marked_state = sparse(marked_state);

	for gam = gammas

		Hb=gam.*Lb-diag(marked_state);
		p_values = zeros(1,length(times));		
		t_count = 1;

		for t = times
			
			p_values(1,t_count) = abs( dot( expm(-i*Hb*t) * initial_state, marked_state ) )^2;
            % careful when using this 'definition' of the optimal gamma (which is actually the
            % correct one) since one has to exclude t = 0 on the one hand but not start to late on
            % the other hand in order to not miss the minimum); also swap comments for the command
            % for max_p_values below
            %p_values(1,t_count) = t/abs( dot( expm(-i*Hb*t) * initial_state, marked_state ) )^2;
			
            t_count = t_count + 1;

		end

		max_p_values( n_count, gam_count ) = max( p_values );
        %%max_p_values( n_count, gam_count ) = min( p_values );
		
        fprintf('%d %f\n', n, gam)
		gam_count = gam_count + 1;

	end

	n_count = n_count + 1;
	t0=toc;
	fprintf('time: %f\n',t0)

end

semilogy(gammas,max_p_values(1,:),'r',gammas,max_p_values(2,:),'m',gammas,max_p_values(3,:),'b')
%semilogy(gammas,max_p_values(1,:),'r')
