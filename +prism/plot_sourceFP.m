function plot_sourceFP(FPs, ROI, FPname)

if nargin <3
    
    FPname = {'FP1', 'FP2'};
end

[nX, nY, nS] = size(FPs); 

dsRedImg = FPs(:,:,1);
mCherryImg = FPs(:,:,2);

if nargin <2 || isempty(ROI)
    ROI = [1 nX 1 nY];
end

s = figure;

rm =  cat(1, [1 1 1], cbrewer('seq', 'Reds',100, 'linear'));
% bm =  cat(1, [1 1 1], cbrewer('seq', 'Blues',100, 'linear'));
% gm =  cat(1, [1 1 1], cbrewer('seq', 'Greens',100, 'linear'));
pm = cat(1, [1 1 1], cbrewer('seq', 'RdPu',100, 'linear'));


lims1 = [0, max(dsRedImg(:))*0.1];
lims2 = [0, max(mCherryImg(:))*0.1];


r1 = subplot(1,3,1);
imagesc(dsRedImg); axis image; caxis(lims1);colorbar;
colormap(r1, rm);
formatAxes
    xlim([ROI(1) ROI(2)]);ylim([ROI(3) ROI(4)]);

title(FPname{1})

r2 = subplot(1,3,2);
imagesc(imgaussfilt(mCherryImg)); axis image;caxis(lims2); colorbar;
colormap(r2, pm);
linkaxes([r1 r2], 'xy')
    xlim([ROI(1) ROI(2)]);ylim([ROI(3) ROI(4)]);

formatAxes
title(FPname{2})


r3 = subplot(1,3,3);

plotOpts.lmax = max(FPs(:));
plotOpts.lmin = min(FPs(:));
plotOpts.bin_size = range([plotOpts.lmin,plotOpts.lmax]/100);
plotOpts.smooth_sigma = 1;

[xy, bins] = xy_projection(reshape(FPs, [], 2), plotOpts);

imagesc(bins, bins, imgaussfilt(xy)); axis image; axis xy
caxis([0 0.00001])
colormap(r3, flip(gray));

xlabel(FPname{1})
ylabel(FPname{2})

xlim([plotOpts.lmin plotOpts.lmax])
ylim([plotOpts.lmin plotOpts.lmax])
formatAxes



end