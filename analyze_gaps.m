%% This script computes the lowest-eigenvalues gap of the full-size
%% Hamiltonian ('f') and of the reduced Hamiltonian ('r') -- and it shows
%% that they are not the same in the root case (which means that in this
%% case the first excited state of the full-size system does not overlap
%% with the intial state and hence the corresponding eigenvalue plays 
%% no part in the dynamics)
%% and that this gap of the reduced Hamiltonian predicts the measurement
%% time correctly in the root case l_mv = 1 (--> predicted time gives
%% success probability p = 0.5 ). 
%% The method fails for cases other than l_mv = 1 (even if one adjusts
%% gamma accordingly) 


% n_max and n_min for looping over different system sizes
n_min = 8;
n_max = 12; % for n_max >= 17 array exceeds maximum array size in Matlab
            % (since we are working with the full operators here)

gamma = 1.0

format short;

% preparations for saving and displaying the data
%n_rows = (n_max+n_min)*(n_max-n_min+1)/2;
n_rows = n_max-n_min+1;
row_counter = 1;
data = zeros(n_rows,10);
fprintf('\n  n  l_mv  ev1_f  ev2_f  t_0_f  succ_prob_f  ev1_r  ev2_r  t_0_r  succ_prob_r  \n \n');

% options for the numerical evaluation of the eigenvalues
opts.isreal = 1;
opts.issym = 1;
opts.tol = 10^(-6);
opts.maxit = 1000;
warning('off', 'MATLAB:eigs:SigmaChangedToSA');

% loop over different system sizes
for num_l = n_min : n_max
    
    % system size and full-size initial state (uniform distribution)
    N=2^num_l-1;
    uni_f = 1 / sqrt( N ) * ones( N, 1);
    
    %for l_mv = 1 : num_l
    for l_mv = 1 : 1
        
        data(row_counter,1)=num_l;
        data(row_counter,2)=l_mv; 
        
        % marked state in the full-size system
        mark_f =  zeros( N, 1);
	    mark_f(2^(l_mv-1),1)=1;

        % generate full Hamiltonian
        H = sparse(gamma.*generate_L(num_l)-diag(mark_f));
        
        % compute eigenvalues
        evs = eigs(H+5*speye(size(H)),[],2,'SR',opts);
        evs = evs-5;
           
        % fill in data
        data(row_counter,3)=min(evs);
        data(row_counter,4)=max(evs);
        data(row_counter,5)=pi/(abs(evs(1)-evs(2)));
        
        % find probability at the time predicted by the gap by evaluating the propagator
        data(row_counter,6)=abs(dot(expm(-i*H*data(row_counter,5)) * uni_f , mark_f))^2;

        % generate reduction matrix and operators/vectors in the reduced system
        V = generate_V(num_l,l_mv);
        Vp = sparse(pinv(V));
        V = sparse(V);
        Hb = V*H*Vp;
        mark_r = V * mark_f;
	    uni_r = V * uni_f;
        
        % compute eigenvalues
        evs = eigs(Hb+5*speye(size(Hb)),[],2,'SR',opts);
	    evs = evs-5;

        % fill in data
        data(row_counter,7)=min(evs);
        data(row_counter,8)=max(evs);
        data(row_counter,9)=pi/(abs(evs(1)-evs(2)));	

        % generate reduction matrix and operators/vectors in the reduced system
        data(row_counter,10)=abs(dot(expm(-i*Hb*data(row_counter,9)) * uni_r , mark_r))^2;
	    
        row_counter = row_counter + 1;

    end

end

disp(data)
