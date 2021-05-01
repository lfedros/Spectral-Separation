function [fpR, fpG] = mixFP(FPs, waveL)
%% INPUTS
%   FPs      1*nFPs cell array containing names of query FPs. Look-up tables for spectral data
%            are stored in the folder 'GitHub/Spectral
%            Separation/2Pspectra'. Source:
%            FPbase.com and Dr. M. Drobizhev

%   waveL    1*nW query excitation wavelenghts
%% OUTPUTS

%   fpR     nW*nFPs vector of mixing coefficients, estimated from in vitro excitation
%           and emission spectra, for  fluorescence recorded in user defined (hard coded) green and
%           red channel
%%

nFP = numel(FPs);

for iF = 1:nFP
    
    fp_name = FPs{iF};
    
    % get excitation spectrum
    excSpectrum(:, iF) = prism.excitationSpectrum(fp_name, waveL);
    
    %get emission spectrum, weighted by bscope red channel 607/70nm
    [emsR(:,iF), ~, pmt] = prism.emissionSpectrum(fp_name, 572:642);

    %get emission spectrum, weighted by bscope green channel 525/50nm
    emsG(:,iF) = prism.emissionSpectrum(fp_name, 500:550);
    
    %calculate mixing coefficient in R and G channel
    fpR(:, iF) =  excSpectrum(:, iF).*nansum(emsR(:,iF));
    fpG(:, iF) =  excSpectrum(:, iF).*nansum(emsG(:,iF));
    
end


end