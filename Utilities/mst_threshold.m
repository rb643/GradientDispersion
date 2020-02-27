function [A] = mst_threshold(Co,MyCost,bin)

% input
% Co = weighted matrix
% MyCost = cost (in range [0,1])
% bin = "binary" = optional flag: T if binary (default), F if weighted
%
% output
% A = thresholded matrix
%
% author: Petra V?rtes
% modified by Frantisek Vasa

if nargin < 3
    bin = true;
end

n=size(Co,1);       % N nodes
Co=(Co+Co')/2;      % force symmetrize matrix
Co(find(Co<0))=0;   % set negative correlations to 0
Co(1:n+1:n*n)=1;    % set diagonal to ones

% create MST (the minimum spanning tree of the network)  
D=ones(n,n)./Co;
MST=kruskal_mst(sparse(D));

% order C according to decreasing wieghts in the correlation matrix
Co=triu(Co,1);
ind = find(triu(ones(n,n),1));
Clist = Co(ind);
Cnonz = length(Clist);
[ClistSort, IX] = sort(Clist,'descend');
[row col]=ind2sub([n,n],ind(IX));
dd = length(Clist);

% store Initial MST in the adjacency matrix A that defines the output network
A = double(logical(full(MST)));

% grow the network according to weights in Co matrix 
t=1;
enum=n-1;
% add edges in correct order until all possible edges exist
while (enum < MyCost*n*(n-1)/2)
    % if edge wasn't initially included in MST
    if A(row(t),col(t)) == 0
        if bin
            %add edge
            A(row(t),col(t)) = 1; %Co(row(t),col(t)); % binary version
            A(col(t),row(t)) = 1; %Co(col(t),row(t)); % binary version
            enum=enum+1;
        else
            A(row(t),col(t)) = Co(row(t),col(t)); % weighted version
            A(col(t),row(t)) = Co(row(t),col(t)); % weighted version
            enum=enum+1;
        end
    end
    t=t+1;
end

end