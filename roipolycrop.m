function [axialROImask, axialROI] = roipolycrop(im, varargin)
% im - input image
% draw circle on image
% return mask
Pos = ParseInputs('Position',[] ,varargin)
[r c wid] = size(im);
if wid==1
    k(:,:,1)=im./max(im(:));
    k(:,:,2)=im./max(im(:));
    k(:,:,3)=im./max(im(:));
elseif wid==3
    k=im./max(im(:));
else
    error('input must be RGB or grayscale')
end
figure(332);
imagesc(k);

axialROI=impoly(gca, Pos);
axialROImask=createMask(axialROI);
return;
