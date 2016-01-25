% This script is similar to 'search_evolution', it generates plots
% of the evolution at the marked vertex. It computes three systems
% at once and each time it gives the plot of p as a function of
% time and the plot of t/p
%
% --> to experiment with the behavior, change the following parameters:
% num_l: size of the system (number of levels)
% l_mv: location of the marked vertex (level on which it is)
% gamma: value of the search parameter (see paper for best values)
% t_max: time up to which the evolution is simulated

clear;

%% FIRST

% set time interval
t_max = 20;
num_t = 200;
times = [0:t_max/num_t:t_max];

% set number of levels and the level of the marked vertex
num_l = 8;
l_mv = 4;

% set gamma (--> use find_gammas.m to find best value) and system size
gamma = 2/3;
N = 2^num_l-1;

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

clf;
subplot(2,3,1)
plot(times,graph_r);
title(sprintf('p for n=%d l=%d (gamma=%f)',num_l,l_mv,gamma))


%% SECOND

% set time interval
%t_max = 500;
%num_t = 200;
times = [0:t_max/num_t:t_max];

% initialize array to store the data
graph_r = zeros(size(times));

% index counter
k=1;

% loop over times and evaluate full and reduced evolution
for t = times
	graph_r(k) = t/abs(dot(expm(-i*Hb*t) * initial_r , marked_r))^2;		
	k = k+1;	
end

subplot(2,3,4)
plot(times,graph_r);
%axis([ 0 t_max 0 250 ])
title(sprintf('t/p for n=%d l=%d (gamma=%f)',num_l,l_mv,gamma))


%% THIRD

% set time interval
t_max = 40;
num_t = 200;
times = [0:t_max/num_t:t_max];

% set number of levels and the level of the marked vertex
num_l = 12;
l_mv = 6;

% set gamma (--> use find_gammas.m to find best value) and system size
gamma = 2/3;
N = 2^num_l-1;

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

subplot(2,3,2)
plot(times,graph_r);
title(sprintf('p for n=%d l=%d (gamma=%f)',num_l,l_mv,gamma))


%% FOURTH

% set time interval
%t_max = 1000;
%num_t = 200;
times = [0:t_max/num_t:t_max];

% initialize array to store the data
graph_r = zeros(size(times));

% index counter
k=1;

% loop over times and evaluate full and reduced evolution
for t = times
	graph_r(k) = t/abs(dot(expm(-i*Hb*t) * initial_r , marked_r))^2;		
	k = k+1;	
end

subplot(2,3,5)
plot(times,graph_r);
%axis([ 0 t_max 0 250 ])
title(sprintf('t/p for n=%d l=%d (gamma=%f)',num_l,l_mv,gamma))


%% FIFTH

% set time interval
t_max = 80;
num_t = 200;
times = [0:t_max/num_t:t_max];

% set number of levels and the level of the marked vertex
num_l = 16;
l_mv = 8;

% set gamma (--> use find_gammas.m to find best value) and system size
gamma = 2/3;
N = 2^num_l-1;

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

subplot(2,3,3)
plot(times,graph_r);
title(sprintf('p for n=%d l=%d (gamma=%f)',num_l,l_mv,gamma))


%% SIXTH

% set time interval
%t_max = 1000;
%num_t = 200;
times = [0:t_max/num_t:t_max];

% initialize array to store the data
graph_r = zeros(size(times));

% index counter
k=1;

% loop over times and evaluate full and reduced evolution
for t = times
	graph_r(k) = t/abs(dot(expm(-i*Hb*t) * initial_r , marked_r))^2;		
	k = k+1;	
end

subplot(2,3,6)
plot(times,graph_r);
%axis([ 0 t_max 0 250 ])
title(sprintf('t/p for n=%d l=%d (gamma=%f)',num_l,l_mv,gamma))
