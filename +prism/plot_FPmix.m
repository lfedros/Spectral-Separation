function plot_FPmix(FPmix,waveL, ROI)

[nX, nY, nW] = size(FPmix);

if nargin <3
    ROI = [1 nX 1 nY];
end

[waveL, sortW] = sort(waveL, 'ascend');

FPmix = FPmix(:,:,sortW);

plotOpts.lmax = max(FPmix(:));
plotOpts.lmin = min(FPmix(:));
plotOpts.bin_size = range([plotOpts.lmin, plotOpts.lmax])/100;

plotOpts.smooth_sigma = 1;

r = figure;

for iW = 1:nW
    
    rimg{iW}= subplot(nW,nW,iW + (iW-1)*(nW));
    imagesc(imgaussfilt(FPmix(:,:,iW))); axis image;
    colormap(rimg{iW}, flip(gray))
    caxis(prctile(makeVec(FPmix(:,:,iW)), [30, 98]));
    
    formatAxes
    set(gca, 'Xtick', [], 'YTick', [])
    xlim([ROI(1) ROI(2)]);ylim([ROI(3) ROI(4)]);
    
    title ([num2str(waveL(iW)), ' nm']);

    for iWW = iW+1:nW
        thisW = FPmix(:,:,iW);
        otherW = FPmix(:,:,iWW);
        [xy, bins] = xy_projection(cat(2, thisW(:), otherW(:)), plotOpts);
        rdist{iW, iWW}= subplot(nW,nW,iWW + (iW-1)*(nW));
        imagesc(bins, bins, imgaussfilt(xy)); axis image; axis xy
        colormap(rdist{iW, iWW}, flip(gray));
        caxis([0 0.00001])

        xlim([plotOpts.lmin plotOpts.lmax])
        ylim([plotOpts.lmin plotOpts.lmax])
        formatAxes
    end
    
    
end



end