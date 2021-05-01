function  [sources, mixing, finalErr] = ...
    learnSources_iterative(data, mixGuess,forceRatio, doPlot)
% learnSources_iterative iterative linear unmixing algorithm to learn both
% the source images and the mixing coefficients, constructed to return maximally 
% uncorrelated source images while minimizing the quadratic reconstruction error over the data. 
% Starting from the theoretical mix- ing coefficients, iterates between estimating the source 
% images while keeping the mixing coefficient fixed and optimizing the mixing coefficients while keeping the source images fixed, until convergence
%%


if nargin < 3
    forceRatio = false;
end

if nargin < 4
    doPlot = false;
    
end
% add intercept?
% mixGuess = cat(2, mixGuess, ones(size(mixGuess, 1),1));

if forceRatio
    r0 = mixGuess(:,2)./mixGuess(:,1);
    r = r0;
end

data = data'; % nPx*nCh --> nCh*nPx

%  estimate sources with initial mixing guess
sourceGuess = pinv(mixGuess)*data; % data is nCh*nPx, mixing is nCh*nF, sources is nF*nPx

% initialise loss function
errReconstruct = data - mixGuess*sourceGuess;
errReconstruct = sum(errReconstruct(:).^2)/sum(data(:).^2); % reconstruction error
uncorr = abs(corr(sourceGuess(1,:)', sourceGuess(2,:)')); % correlation between the 2 sources

% negWeigth = sum(sourceGuess(sourceGuess(2,:)<0).^2)...
%     /sum(sourceGuess(2, :).^2); % penalty for non-sensical negative values in source
% mutInfo = mutualInformation(sourceGuess(1,:), sourceGuess(2,:)); % mutual information between the sources
% sourceBlow = sum(sourceGuess(2,:).^2)/sum(sourceGuess(:).^2); % penalty to avoid one source to become dominant;

iter = 1;
if forceRatio
    err(iter)  = errReconstruct*uncorr;
else
    err(iter)  = errReconstruct*uncorr;
    
end
tol = err(iter);

if doPlot
    e = figure;
end

finalErr = nan(100,1);
while  tol > 10^(-10) && iter <100
    
    iter = iter + 1;
    
    % improve guess of mixing matrix
    mixing = pinv(sourceGuess')*data'; % data is nPx *nCh, sources is nPx*nF, mixing is nF*nCh
    mixing = mixing';
    mixing(mixing<0) = 0.000001; % negative weights do not make sense
    
    if forceRatio
        % allow mixing coefficient to scale, but constrain relative ratio to be within range from initial guess
        r = mixing(:,2)./mixing(:,1);
        for iW = 1:size(mixing,1)
            if r0(iW)/r(iW) <0.1 || r0(iW)/r(iW) >10 % if relative ratio changes too much, reset
                %             mixing(iW, 2:end) =  mixing(iW, 1)*mixRatio(iW, 2:end);
                mixing(iW, 2) =  mixing(iW, 1)*r0(iW);
            end
            %        mixing(:, 2) = bsxfun(@times, mixing(:,1),mixRatio(:,2));
        end
    end
    % recompute sources
    sources = pinv(mixing)*data;
    
    % calculate loss function
    
    errReconstruct = data - mixing*sources;
    errReconstruct = sum(errReconstruct(:).^2)/sum(data(:).^2);
    uncorr = abs(corr(sourceGuess(1,:)', sourceGuess(2,:)'));
    %     negWeigth = sum(sources(sources(2,:)<0).^2)/sum(sources(2, :).^2);
    %     mutInfo = mutualInformation(sources(1,:), sources(2,:));
    %     sourceBlow = sum(sources(2,:).^2)/sum(sources(:).^2);
    
    if forceRatio
        err(iter)  = errReconstruct*uncorr;
    else
        err(iter)  = errReconstruct*uncorr;
        
    end
    
    % compute improvement from previous iteration
    tol =  err(iter-1) - err(iter);
    
    % plot loss function if requested
    if doPlot
        figure(e);
        cla;plot(1:iter, log10(err), '-o'); hold on; drawnow;
        xlabel('iteration #')
        ylabel('log10(err)')
    end
    % pause;
    
    % break iterations when minimum is found or tolerance criterion is met
    if tol < 0 % if minimum was previous iteration, go back one step
        mixing = mixGuess;
        sources = sourceGuess;
        finalErr(1:iter-1) = err(1:iter-1);
    else % if there was improvement, record it
        mixGuess = mixing;
        sourceGuess = sources;
        finalErr(1:iter) = err(1:iter);
    end
    
end

sources = pinv(mixing)*data;



end