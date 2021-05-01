
clear;

% set paths to code
prism_repo = '/Users/lfedros/Documents/GitHub/Spectral Separation';
addpath(genpath(prism_repo));

fedbox_repo = '/Users/lfedros/Documents/GitHub/FedBox';
addpath(genpath(fedbox_repo));

%% create test source images

source_dsRed = rgb2gray(insertText(zeros(512,512), [5,100], 'dsRed', ...
    'TextColor', 'white', 'BoxColor', 'black', 'FontSize', 150));

source_mCherry = rgb2gray(insertText(zeros(512,512), [5,150], 'mCherry', ...
    'TextColor', 'white', 'BoxColor', 'black', 'FontSize', 110));

%% mix them

% these are the wavelengths we imaged
waveL = [780 890 970 1020];

% calculate mixing coefficient according to excitation and emission spectra
[fpR, fpG] = prism.mixFP({'dsRed2', 'mCherry'}, waveL);

nW = numel(waveL);
[nX, nY] = size(source_dsRed);

% add some noise to the mixing coefficient
noisy_fpR = fpR + randn(size(fpR)).*fpR/50;

% create synthetic spectral data
for iW = 1:nW
    % add noise to each image
    noise = randn(512)*max(fpR(iW,:))/50;
    mixed_red(:,:,iW) = noisy_fpR(iW,1).* source_dsRed + noisy_fpR(iW,2).* source_mCherry + noise;
    
end

% plot the data
prism.plot_FPmix(mixed_red, waveL);

%% spectral demixing

% we are going to use all the pixels in the image
forePx = reshape(mixed_red, nX*nY, nW);

% demix
[sources, mixing, iterErr] = prism.learnSources_iterative(forePx, fpR,0, 1);

% prepare source images for plotting
dsRedImg = reshape(sources(1,:), nX, nY);
mCherryImg= reshape(sources(2,:), nX, nY);

% plot
prism.plot_sourceFP(cat(3, dsRedImg,  mCherryImg));

%% now let's try with some real data
clear;

% set paths to code
repo_prism = '/Users/lfedros/Documents/GitHub/Spectral Separation';
addpath(genpath(repo_prism));

fedbox = '/Users/lfedros/Documents/GitHub/FedBox';
addpath(genpath('/Users/lfedros/Documents/GitHub/FedBox'));

% load the data
load('FR140_exp03_plane04');

% estimate mixing coefficient based on excitation and emission spectra
[fpR, fpG] = prism.mixFP(FPs, waveL);

[nX, nY, nW] = size(mixed_red);

% plot the data
prism.plot_FPmix(mixed_red, waveL, [170, 470, 80, 380]);

% segment the foreground neurons (labelling is very sparse, we want to use
% only pixels with some signal)

bw = prism.foreground_bw(mixed_red);

forePx = reshape(mixed_red, nX*nY, nW);

forePx = forePx(bw(:), :);

% iterative demixing
[~, mixing, iterErr] = prism.learnSources_iterative(forePx, fpR,1, 1);
sources = pinv(mixing)*reshape(mixed_red, nX*nY, nW)';

% plot the results
dsRedImg = reshape(sources(1,:), nX, nY);
mCherryImg= reshape(sources(2,:), nX, nY);

prism.plot_sourceFP( cat(3, dsRedImg,  mCherryImg), [170, 470, 80, 380]);

