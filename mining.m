%%
%  Clustering gene expression data
%  `https://www.mathworks.com/help/bioinfo/examples/`
%%


cd('~/Desktop');
f = fopen('fit.csv');
df = textscan(f, '%s %f %f %f', 'Delimiter', ',');
fclose(f);

genes = df{1};
exprs = horzcat(df{2:4});

%% filter gene variances
mask = genevarfilter(exprs, 'Percentile', 60);
genes = genes(mask);
exprs = exprs(mask,:);

%% calculate distances based on correlation
corrDist = pdist(exprs,'corr');
clusterTree = linkage(corrDist,'average');
clusters = cluster(clusterTree,'maxclust',8);

%% plot
figure
for c = 1:16
    subplot(4,4,c);
    plot(exprs((clusters == c),:)');
    axis tight
end
suptitle('Hierarchical Clustering of Profiles');

%% k-means clustering
rng('default');

[cidx, ctrs] = kmeans(exprs,9,'dist','corr','rep',5,'disp','final');
figure
for c = 1:16
    subplot(4,4,c);
    plot(exprs((cidx == c),:)');
    axis tight
end
suptitle('K-Means Clustering of Profiles');

%% plot
figure
for c = 1:16
    subplot(4,4,c);
    plot(ctrs(c,:)');
    axis tight
    axis off
end
suptitle('K-Means Clustering of Profiles');

%% plot PCA
[pc, zscores, pcvars] = pca(exprs);
figure
scatter(zscores(:,1),zscores(:,2));
xlabel('PC1');
ylabel('PC2');
title('Principal Component Scatter Plot');
figure

pcclusters = clusterdata(zscores(:,1:2),'maxclust',7,'linkage','av');
gscatter(zscores(:,1),zscores(:,2),pcclusters)
xlabel('First Principal Component');
ylabel('Second Principal Component');
title('Principal Component Scatter Plot with Colored Clusters');

B = unique(pcclusters);
out = [B,histc(pcclusters,B)]
