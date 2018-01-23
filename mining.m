cd('~/Desktop');
f = fopen('fit2.csv');
df = textscan(f, '%s %f %f %f', 'Delimiter', ',');
fclose(f);

genes = df{1};
exprs = horzcat(df{2:4});

%% gene var filter
mask = genevarfilter(exprs, 'Percentile', 60);
genes = genes(mask);
exprs = exprs(mask,:);

%% 
corrDist = pdist(exprs,'corr');
clusterTree = linkage(corrDist,'average');
clusters = cluster(clusterTree,'maxclust',8);

%%
figure
for c = 1:16
    subplot(4,4,c);
    plot(exprs((clusters == c),:)');
    axis tight
end
suptitle('Hierarchical Clustering of Profiles');

%%
rng('default');

[cidx, ctrs] = kmeans(exprs,9,'dist','corr','rep',5,'disp','final');
figure
for c = 1:16
    subplot(4,4,c);
    plot(exprs((cidx == c),:)');
    axis tight
end
suptitle('K-Means Clustering of Profiles');

%%
figure
for c = 1:16
    subplot(4,4,c);
    plot(ctrs(c,:)');
    axis tight
    axis off
end
suptitle('K-Means Clustering of Profiles');

%% 
[pc, zscores, pcvars] = pca(exprs);
figure
scatter(zscores(:,1),zscores(:,2));
xlabel('First Principal Component');
ylabel('Second Principal Component');
title('Principal Component Scatter Plot');
figure

%%
pcclusters = clusterdata(zscores(:,1:2),'maxclust',7,'linkage','av');
gscatter(zscores(:,1),zscores(:,2),pcclusters)
xlabel('First Principal Component');
ylabel('Second Principal Component');
title('Principal Component Scatter Plot with Colored Clusters');
%%
B = unique(pcclusters);
out = [B,histc(pcclusters,B)]

%% 
mask = pcclusters == 3;
genes = genes(mask);
exprs = exprs(mask,:);





