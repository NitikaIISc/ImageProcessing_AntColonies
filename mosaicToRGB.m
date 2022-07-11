function img_out = mosaicToRGB(im)    
    G = im(1:2:end,2:2:end);
    R = im(1:2:end,1:2:end);
    B = im(2:2:end,2:2:end);
    img_out = cat(3,R,G,B);
end