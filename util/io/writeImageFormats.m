function writeImageFormats(data, spacing, window, fileBaseName)

% store as MHD+RAW
writeMetaImage(data, [fileBaseName '.mhd'], spacing);

% store as PNG
hWindowFunction = @(val) (val-window(1))/(window(2)-window(1));
imwrite(hWindowFunction(data'), [fileBaseName '.png']);
