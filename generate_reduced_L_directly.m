%% generates the graph Laplacian on the reduced system
%% (--> comb-structure) directly (without generating the full
%% graph Laplacian and the reduction matrix first)
%
% input:
% num_levels --> number of levels of the tree that
%                is to be reduced
% level_mv --> level on which the marked vertex (mv) is
%              (we always let the marked vertex be the
%               the leftmost node on the given level)
%
% output:
% Lb --> the graph Laplacian on the comb (--> "L bar")
% imv --> index (in the comb) of the marked vertex
% marked_state --> vector |w> in the reduced system
% initial_state --> vector |s> in the reduced system


function [ Lb, imv, marked_state, initial_state ] = generate_reduced_L_directly( num_levels, level_mv )

	function biggest_child = biggestChild(index)
		indexes = find(Lb(index,:));
		biggest_child = indexes(end);
		if biggest_child <= index
			biggest_child = 0;
		end
	end

	function par = parent(index)
		if index == 1
			par = 0;
		else		
			indexes = find(Lb(index,:));
			par = indexes(1);
		end
	end


% preparations
N = 2^num_levels-1;                                   % size of the orig system
dim = - level_mv^2/2 + level_mv*( num_levels + 1/2 ); % size of the comb
imv = level_mv^2/2 - level_mv/2 + 1;                  % index of the marked vertex
Lb = zeros( dim, dim );                               % intialize Lb


% Lb: diagonal entries
for k = 1 : dim
	Lb( k, k ) = 3;
end
for k = 1 : level_mv
	Lb( dim - k + 1, dim - k + 1 ) = 1;
end 
Lb(1,1) = 2;


% Lb: other entries in the root case
if level_mv == 1
	offdiag = sqrt(2).*ones(dim-1,1);
	Lb = Lb-diag(offdiag,+1)-diag(offdiag,-1);

% Lb: other entries in general
else
	for k = dim:-1:imv+level_mv
		Lb(k,k-level_mv)=-sqrt(2);
		Lb(k-level_mv,k)=-sqrt(2);	
	end
	next = imv+level_mv-1;
	for r = level_mv-1:-1:1
		for k = next:-1:next-r+2
			Lb(k,k-r-1)=-sqrt(2);
			Lb(k-r-1,k)=-sqrt(2);	
		end
		Lb(next-r+1,next-2*r)=-1;
		Lb(next-2*r,next-r+1)=-1;	
		Lb(next-r,next-2*r)=-1;
		Lb(next-2*r,next-r)=-1;
		next = next-r-1;
	end	

end


% marked_state
marked_state = zeros(dim,1);
marked_state( imv,1) = 1;


% initial_state in the root case and in general
initial_state = ones(dim,1);
if level_mv == 1
	for k = 2:dim
		initial_state(k,1)=2^(k-1);
	end
else
	curr = biggestChild(imv);
	p = 1;
	while  curr~=0
		initial_state(curr)=2^p;
		p = p+1;
		curr = biggestChild(curr);
	end
	par = parent(imv);
	while par~=0
		curr = biggestChild(par);
		p = 0;
		while  curr~=0
			initial_state(curr)=2^p;
			p = p+1;
			curr = biggestChild(curr);
		end
		par = parent(par);
	end
end
initial_state = 1/sqrt(N).*sqrt(initial_state);		


end
