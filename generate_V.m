%% generates the reduction matrix V by grouping certain
%% nodes of the original (full-size) tree together
%
% input:
% num_levels --> number of levels of the tree that
%                is to be reduced
% level_mv --> level on which the marked vertex (mv) is
%              (we always let the marked vertex be the
%               the leftmost node on the given level)
%
% output:
% V --> the reduction matrix


function [ V ] = generate_V( num_levels, level_mv )


%% auxiliary functions

%% appends a cell to an array of cells
%
% (used by 'generate_V'; in this context the cells are indices
% of the nodes that will be grouped together by the reduction)
function [ cellArray ] = appendTo( cellArray, element )
	cellArray{ length( cellArray ) + 1 } = element;
end

%% sorts an array of lists of numbers (cell array) according
%% to the first entry of each list/cell
%
% (used by 'generate_V'; in this context the cells are indices
% of the nodes that will be grouped together by the reduction)
function [ sortedArray ] = sortCellArray( cellArray )
	d = length(cellArray);
	temp = zeros(1,d);
	for k = 1:d
		temp(k) = cellArray{k}(1 );
	end
	[ordered,permutation] = sort(temp);
	sortedArray = {};
	for k = 1 : length( permutation )
		sortedArray = appendTo( sortedArray, cellArray{ permutation( k ) } );
	end	
end


%% code for generate_V


N = 2^num_levels - 1;	    % total number of nodes
index_mv = 2^(level_mv-1);  % index of the mv
groups = {};         	    % pool of groups
temp_groups = {};    	    % temporary pool of groups
% temp_groups always contains the nodes from which we
% will continue in the next step

% add all one-element groups to pool of groups
% --> all ancestors of the mv ...
for k = 1 : level_mv - 1
    groups = appendTo( groups, 2^(k-1) );
    groups = appendTo( groups, 2^(k)+1 );
    temp_groups = appendTo( temp_groups, 2^(k)+1 );
end
% ... and the mv itself
groups = appendTo( groups, index_mv );
temp_groups = appendTo( temp_groups, index_mv );

% add groups with j^2 elements to the pool of groups
for j = 1 : num_levels
    curr_groups = temp_groups;
    temp_groups = {};    
    for curr = curr_groups
        next_elem = 2*curr{1};
        if next_elem<=N
            next_group = next_elem : next_elem+2^j - 1;
            groups = appendTo( groups, next_group );
            temp_groups = appendTo( temp_groups, next_group);
	end
    end
end

% now generate the reduction matrix V, which is of 
% size length(groups) x N    
V = zeros( length( groups ), N );
groups = sortCellArray(groups);
columnCounter = 1;
for row = 1 : length( groups )
	num_curr_group = length( groups{ row } );
	for column = 1 : num_curr_group
		V( row, columnCounter ) = 1 / sqrt( num_curr_group );
		columnCounter = columnCounter + 1;		
	end
end

end
