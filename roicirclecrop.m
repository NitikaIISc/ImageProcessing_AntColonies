function [axialROImask, axialROI] = roicirclecrop(im, varargin)
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

axialROI=imellipse(gca, Pos);
%axialROI=impoly(gca, Pos);
pause
axialROImask=createMask(axialROI);
return;
