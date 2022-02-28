%% DEMO FOR ITERATIVE, LINEAR SPECTRAL UNMIXING
% synthetic data

%% DEMO FOR ITERATIVE, LINEAR SPECTRAL UNMIXING
% Data from Rossi et al. 2020, Nature. 

clear;

% set paths to code
repo_prism = '/Users/lfedros/Documents/GitHub/Spectral Separation';
addpath(genpath(repo_prism));

fedbox = '/Users/lfedros/Documents/GitHub/FedBox';
addpath(genpath('/Users/lfedros/Documents/GitHub/FedBox'));
addpath('/Users/lfedros/Documents/GitHub/cameraLucida/Utils');

% load the data
load('/Users/lfedros/Downloads/test_SB_790_960_gcamp.mat');
load('/Users/lfedros/Downloads/test_SB_790_960.mat');

%% subtract gcamp

R = mixed_red;
G = mixed_green;

 [curedR, bleedThrough, redList, redAll, globalThresh] = ...
    bleedCureNL(R, G, [], 5, [], 1);

mixed_red = curedR;
FPs = {'tdTomato', 'mCherry'};
waveL = [790, 960];
%% estimate mixing coefficient based on excitation and emission spectra
[fpR, fpG] = prism.mixFP(FPs, waveL);

[nX, nY, nW] = size(mixed_red);

%% plot the data
prism.plot_FPmix(mixed_red, waveL);

%% segment the foreground neurons 
% labelling is very sparse, we want to use only pixels with some signal)

bw = prism.foreground_bw(mixed_red);

forePx = reshape(mixed_red, nX*nY, nW);

forePx = forePx(bw(:), :);

%% iterative, spectral linear unmixing
[~, mixing, iterErr] = prism.learnSources_iterative(forePx, fpR,0, 1);
sources = pinv(mixing)*reshape(mixed_red, nX*nY, nW)';

dsRedImg = reshape(sources(1,:), nX, nY);
mCherryImg= reshape(sources(2,:), nX, nY);

%% plot the results

prism.plot_sourceFP( cat(3, dsRedImg,  mCherryImg), [], FPs);

