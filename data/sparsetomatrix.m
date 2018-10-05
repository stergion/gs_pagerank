function [ matrix ] = sparsetomatrix( dirname )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
    
    nodesfile = sprintf('%s%s',dirname,'/nodes');
    fnodes = fopen(nodesfile, 'r');
    Nnodes = str2double(fgets(fnodes));
    fclose(fnodes);
    
    sparsefile = sprintf('%s%s',dirname,'/adj_sparse');
    fsparse = fopen(sparsefile, 'r');
    
    sparse = fscanf(fsparse, "%d %d\n");
    matrix = zeros(Nnodes,Nnodes);
    
    for i= 1:2:length(sparse)
        matrix(sparse(i)+1, sparse(i+1)+1) = 1;
    end
    fclose(fsparse);

end

