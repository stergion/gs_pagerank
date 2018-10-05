function [ pr ] = getPageRanks( matrix, savedir )
%GETPAGERANKS Summary of this function goes here
%   Detailed explanation goes here

% Calculate PageRanks
G = digraph(matrix);
pr = centrality(G, 'pagerank', 'MaxIterations', 200, 'FollowProbability', 0.85);

% Save results to text file
prfile = sprintf('%s%s',savedir,'/matlab_pr');
fpr = fopen(prfile, 'w');
fprintf(fpr,'%f\n',pr);
fclose(fpr);

end
