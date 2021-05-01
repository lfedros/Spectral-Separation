
function [bw, idx, w_tresh] = foreground_bw(mixing_red)

[nX, nY, nW] = size(mixing_red);

bw = zeros(nX, nY);
w_tresh =nan(1, nW);

for iW= 1:nW
    
    w_img = mixing_red(:,:, iW); 
    
    w_tresh(iW) = std(w_img(:)); 
    
    bw = w_img > w_tresh(iW) | bw; % & w_img <3000 ; 
    
end

idx = find(bw);


end