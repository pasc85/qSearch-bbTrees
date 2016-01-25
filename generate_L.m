%% generates the full-size graph Laplacian L
%
% input:
% num_levels --> number of levels of the tree that
%                is to be reduced
%
% output:
% L --> the graph Laplacian on the tree


function [ L ] = generate_L( num_levels )

	N = 2^num_levels - 1;

	L = zeros( N );

	L(1,1) = 2;
	for k=2:2:N-1
	    L(k/2,k)=-1;
	    L(k/2,k+1)=-1;
	    L(k,k/2)=-1;
	    L(k+1,k/2)=-1;
	end
	for k=2:(N-1)/2
	    L(k,k)=3;
	end
	for k=(N-1)/2+1:N
	    L(k,k)=1;
	end


end
