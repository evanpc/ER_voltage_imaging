


function tif=readinim(impath)


%loading voltage imaging tif
tif.info=imfinfo(impath);%for recording details about image
tif.xyz=[tif.info(1).Height,tif.info(1).Width,...
    size(tif.info,1)];%number of pix in x/y, and number of z plains (used later)


%making matrix from entire tif and reading images
tif.reads=zeros(tif.xyz(1),tif.xyz(2),tif.xyz(3),'uint16'); %this is the matrix for all im data.

for i=1:tif.xyz(3)
    tif.reads(:,:,i)=imread(fullfile(impath),i);
end


%making ref image
tif.ref=zeros(tif.xyz(1),tif.xyz(2),'uint16');
tif.ref=imadjust(max(tif.reads,[],3));

end