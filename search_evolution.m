% flag --> set to '1' to compare evolution of the full-size and 
% the reduced system for small size (<= 10), and to any other number to 
% evolve large systems using the improved iterative methods
flag = 1;

% set time interval
t_max = 200;
num_t = 300;
times = [0:t_max/num_t:t_max];

% set number of levels and the level of the marked vertex
num_l = 8;
l_mv = 2;

% set gamma (--> use find_gammas.m to find best value) and system size
gamma = 1.0;
N = 2^num_l-1;

if flag == 1

    % generate full operator and vectors
    marked_f =  zeros( N, 1);
    marked_f(2^(l_mv-1),1)=1;
    marked_f=sparse(marked_f);
    initial_f = sparse(1 / sqrt( N ) * ones( N, 1));
    H = sparse(gamma.*generate_L(num_l)-diag(marked_f));

    % generate reduced operator and vectors directly
    [ Lb, imv, marked_r, initial_r ] = generate_reduced_L_directly( num_l, l_mv );
    marked_r=sparse(marked_r);
    initial_r=sparse(initial_r);
    Hb = sparse(gamma.*Lb-diag(marked_r));

    % initialize arrays to store the data
    graph_f = zeros(size(times));
    graph_r = zeros(size(times));

    % index counter
    k=1;

    % loop over times and evaluate full and reduced evolution
    for t = times

	    graph_f(k) = abs(dot(expm(-i*H*t) * initial_f , marked_f))^2;
	    graph_r(k) = abs(dot(expm(-i*Hb*t) * initial_r , marked_r))^2;		
	    k = k+1;	

    end

    % plots
    subplot(1,2,1);
    plot(times,graph_f);
    subplot(1,2,2);
    plot(times,graph_r);

else

    % generate reduced operator and vectors directly
    [ Lb, imv, marked_r, initial_r ] = generate_reduced_L_directly( num_l, l_mv );
    marked_r=sparse(marked_r);
    initial_r=sparse(initial_r);
    Hb = sparse(gamma.*Lb-diag(marked_r));

    % initialize array to store the data
    graph_r = zeros(size(times));

    % index counter
    k=1;

    % loop over times and evaluate full and reduced evolution
    for t = times

	    graph_r(k) = abs(dot(expm(-i*Hb*t) * initial_r , marked_r))^2;		
	    k = k+1;	

    end

    % plots
    clf;
    plot(times,graph_r);

end
