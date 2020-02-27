function A_rand=randomize_matrix(A)

%This code creates a random undirected network from the connectivity
%distribution of an undirected adjacency matrix, ie, the intital matrix
%must be symmetric.

% INPUTS:
%   A: an undirected adjacency matrix (symmetric) with no self connections

% OUTPUTS:
%   A_rand: a comparable random network with same number of nodes and
%       connectivity distribution

% written by Sarah F. Muldoon

num_nodes=length(A);
A_rand=zeros(num_nodes);
mask=triu(ones(num_nodes),1);
grab_indices=find(mask > 0);

orig_edges=A(grab_indices);
num_edges=length(orig_edges);

rand_index=randperm(num_edges);
randomized_edges=orig_edges(rand_index);

edge=1;
    for i=1:num_nodes-1
        for j=i+1:num_nodes
            A_rand(i,j)=randomized_edges(edge);
            A_rand(j,i)=randomized_edges(edge);
            edge=edge+1;
        end
    end
end