%% general notes %%




%FUNCTION:

%this script imports a single color multiplane .tif image, extracts the 
%average fluorescence from a hand drawn ROI over each .tif plane (i.e. time), 
%subtracts the average fluorescence from a hand drawn background ROI, and 
% plots the net fluorescence vs time.

%this script can be done iteratively without reloading

%this script creates a variables and figures directory in the home 
%directory of the .tif image being analyzed. In these directories the 
%workspace and fluorescence vs time image are saved, respectively. the 
%output variable for background subtracted fluorescence is f.fluo




%DATA REQUIREMENTS:

%a .tiff image is all



%INSTRUCTIONS:

%in general, follow instructions on UI boxes and check command line feedback. More explanation below:

%select directory where .tif image is stored
%draw cell ROI as many times as needed until sufficient 
%do the same with background ROI
%change input parameters on GUI as needed
%choose to repeat this process with another file or end




%% loading all files and making directories %%


disp(' ');
disp(' ');
disp(' ');
disp('select recording directory');
disp(' ');
disp(' ');
disp(' ');


%getting recording directory to work within
path.exptdir = uigetdir;
[path.filepath,path.name,~] = fileparts(path.exptdir);
path.expt = fullfile(path.filepath,path.name);


%loading tif
tif.name=dir(fullfile(path.expt,'*.tif')); %finding tif in directroy
tif.info=imfinfo(fullfile(tif.name.folder,tif.name.name));%for recording details about image
tif.xyz=[tif.info(1).Height,tif.info(1).Width,size(tif.info,1)];%number of pix in x/y, and number of z plains (used later)


%making matrix from entire tif and reading images
tif.reads=zeros(tif.xyz(1),tif.xyz(2),tif.xyz(3),'uint16'); %this is the matrix for all im data.

for i=1:tif.xyz(3)
    tif.reads(:,:,i)=imread(fullfile(tif.name.folder,tif.name.name),i);
end


%making ref image
tif.ref=zeros(tif.xyz(1),tif.xyz(2),'uint16');
tif.ref=imadjust(max(tif.reads,[],3));


%command line feedback
%disp(strcat('phys file is:    "',v.name.name,'"'));
disp(strcat('voltage imaging file selected is:    "',tif.name.name,'"'));
disp(' ');
disp(' ');
disp(' ');


%making directories and paths for outputs
mkdir(fullfile(path.expt,'figures'));
path.figs=fullfile(path.expt,'figures');
mkdir(fullfile(path.expt,'variables'));
path.vars=fullfile(path.expt,'variables');



%% making masks for analysis and double checking them %%


%making roi mask
choose=-1;
figs.tit=sprintf('draw cell roi');
figure('Color','white')
while choose<0
    imshow(tif.ref);
    title(figs.tit, 'FontSize', 14);
    roi.frm = imfreehand;
    setColor(roi.frm,'r');
    
    choose = questdlg('cell roi correct?  ', 'check cell roi!','yes','no','yes');
    switch choose
        case 'yes'
            disp('   cell roi set')
            choose = 1;
        case 'no'
            figs.tit=sprintf(' redraw cell roi');
            choose =-1;
    end
end
roi.mask = createMask(roi.frm);%this var may need to not exist
savefig(fullfile(path.figs,strcat(path.name,'_cell_roi.fig')))


%making backgound mask
figs.tit=sprintf('draw background roi');
while choose<2
    imshow(tif.ref);
    title(figs.tit, 'FontSize', 14)
    bk.frm = imfreehand;
    setColor(bk.frm,'c');
    
    choose = questdlg('background roi correct?  ',...
        'check background roi!','yes','no','yes');
    switch choose
        case 'yes'
            choose = 2;
            disp('   background set')
        case 'no'
            figs.tit=('redraw background roi');
            choose =-1;
    end
end
bk.mask = createMask(bk.frm);%this var may need to not exist
savefig(fullfile(path.figs,strcat(path.name,'_bkgrnd_roi.fig')))



%% input parameters %%


figs.dlg.tit='input image sample period';
    figs.dlg.prompt={'image sample period (ms)'};
    figs.dlg.title=figs.dlg.tit;
    figs.dlg.dims=[1 70];
    figs.dlg.definput={'100'};
    paras=str2double(inputdlg(figs.dlg.prompt,figs.dlg.title,...
        figs.dlg.dims,figs.dlg.definput));



%% reading image and pairing data with phys%%


%reading rois on tif.reads matrix
tif.temp=zeros(tif.xyz(1),tif.xyz(2));
roi.fluo=zeros(1,tif.xyz(3));
bk.fluo=zeros(1,tif.xyz(3));
for i=1:tif.xyz(3)
    tif.temp=tif.reads(:,:,i);
    roi.fluo(i)=mean(mean(tif.temp(roi.mask))); %mean of roi on image
    bk.fluo(i)=mean(mean(tif.temp(bk.mask)));%mean of background on image
end
f.fluo=roi.fluo-bk.fluo;%this is fluo for experiment but not compatible with phys sample rate


%set axis scale and get conversion for pairing
tscale = questdlg('select time units for x axis on figures  ',...
    'timescale for figures and variables!','ms','sec','min','sec');
switch tscale
    case 'ms'
        figs.xtitle='miliseconds';
        tscale=paras(1);
    case 'sec'
        figs.xtitle='seconds';
        tscale=paras(1)/1000;
    case 'min'
        figs.xtitle='minutes';
        tscale=paras(1)/60000;
end




%% outputs




%figs
plot(f.fluo(1, :), 'Color',[0.9 0 0])
box off
ax=gca;
ax.FontSize = 14;
xlabel(figs.xtitle, 'FontSize',18)
ylabel('fluorescence', 'FontSize',18)
savefig(fullfile(path.figs,strcat(path.name,'_fluorescence.fig')))

%saving workspace%
save(fullfile(path.vars,strcat(path.name,'_workspace.mat')));



%% end of one recording trial analysis %%
%from here on one cell is complete and all files are saved. However, if other
%recordings may share parameters or masks, the following loop will be
%reusable as many times as needed. In addition, all copying and pasting the
%following sections into old completed workspaces can be done to reuse old
%rois or parameters.

%% Loop for new cells if needed %%

choose = questdlg('process annother recording trial with these parameters?  ',...
    'continue or complete','yes','no','no');
switch choose
    case 'yes'
        disp(' ');
        disp(' ');
        disp(' ');
        disp('select next recording directory');
        disp(' ');
        disp(' ');
        disp(' ');
        choose=-1;
    case 'no'
        disp(' analysis complete');
        choose=100;
end



%% master loop for entire redo %%


while choose<100
    %% loading all files and making directories %%
    
    
    %getting recording directory to work within
    path.exptdir = uigetdir(path.filepath);
    [path.filepath,path.name,~] = fileparts(path.exptdir);
    path.expt = fullfile(path.filepath,path.name);
    
    
    %loading tif info
    tif.name=dir(fullfile(path.expt,'*.tif')); %finding tif in directroy
    tif.info=imfinfo(fullfile(tif.name.folder,tif.name.name));%for recording details about image
    tif.xyz=[tif.info(1).Height,tif.info(1).Width,size(tif.info,1)];
    
    
    %making matrix from entire tif and reading images
    tif.reads=zeros(tif.xyz(1),tif.xyz(2),tif.xyz(3),'uint16'); %this is the matrix for all im data.
    
    for i=1:tif.xyz(3)
        tif.reads(:,:,i)=imread(fullfile(tif.name.folder,tif.name.name),i);
    end
    
    
    %making ref image
    tif.ref=imadjust(max(tif.reads,[],3));
    
    
    %command line feedback

    disp(strcat('voltage imaging file selected is:    "',tif.name.name,'"'));
    disp(' ');
    disp(' ');
    disp(' ');
    
    
    %making directories and paths for outputs
    mkdir(fullfile(path.expt,'figures'));
    path.figs=fullfile(path.expt,'figures');
    mkdir(fullfile(path.expt,'variables'));
    path.vars=fullfile(path.expt,'variables');
    
    
    %% double checking masks %%
    
    
    %check roi mask on new ref im
    figure('Color','white')
    imshow(tif.ref);
    title('logged cell roi on new cell recording', 'FontSize', 14)
    hold on
    visboundaries(roi.mask,'Color','r','LineWidth', 0.5);
    hold off
    
    choose = questdlg('reuse cell roi?  ', 'check cell roi!','yes','no','yes');
    switch choose
        case 'yes'
            disp('   cell roi set')
            savefig(fullfile(path.figs,strcat(path.name,'cell_roi.fig')))
            choose=1;
        case 'no'
            figs.tit=sprintf(' draw new cell roi');
            choose =-1;
            while choose<0
                imshow(tif.ref);
                title(figs.tit, 'FontSize', 14);
                roi.frm = imfreehand;
                setColor(roi.frm,'r');
                
                choose = questdlg('cell roi correct?  ',...
                    'check cell roi!','yes','no','yes');
                switch choose
                    case 'yes'
                        disp('   cell roi set')
                        choose = 1;
                    case 'no'
                        figs.tit=sprintf(' redraw cell roi');
                        choose =-1;
                end
            end
            roi.mask = createMask(roi.frm);%this var may need to not exist
            savefig(fullfile(path.figs,strcat(path.name,'_cell_roi.fig')))
    end
    
    %check background mask on new ref im
    imshow(tif.ref);
    title('logged background roi on new cell recording', 'FontSize', 14)
    hold on
    visboundaries(bk.mask,'Color','c','LineWidth', 0.5);
    hold off
    
    
    choose = questdlg('reuse background roi?  ',...
        'check background roi!','yes','no','yes');
    switch choose
        case 'yes'
            disp('   background roi set')
            choose = 1;
            savefig(fullfile(path.figs,strcat(path.name,'bkgrnd_roi.fig')))
        case 'no'
            figs.tit=sprintf('draw new background roi');
            choose =-1;
            while choose<0
                imshow(tif.ref);
                title(figs.tit, 'FontSize', 14)
                bk.frm = imfreehand;
                setColor(bk.frm,'c');
                
                choose = questdlg('background roi correct?  ',...
                    'check background roi!','yes','no','yes');
                switch choose
                    case 'yes'
                        choose = 1;
                        disp('   background set')
                    case 'no'
                        figs.tit=('redraw background roi');
                        choose =-1;
                end
            end
            bk.mask = createMask(bk.frm);%this var may need to not exist
            savefig(fullfile(path.figs,strcat(path.name,'_bkgrnd_roi.fig')))
            
    end
    
    
    
    
    %% input parameters %%
    
    
   
    
    %% reading image and pairing data with phys%%
    
    
    %reading rois on tif.reads matrix
    tif.temp=zeros(tif.xyz(1),tif.xyz(2));
    roi.fluo=zeros(1,tif.xyz(3));
    bk.fluo=zeros(1,tif.xyz(3));
    for i=1:tif.xyz(3)
        tif.temp=tif.reads(:,:,i);
        roi.fluo(i) = mean(mean(tif.temp(roi.mask))); %mean of roi on image
        bk.fluo(i) = mean(mean(tif.temp(bk.mask)));%mean of background on image
    end
    f.fluo=roi.fluo-bk.fluo;%this is fluo for experiment but not compatible withh phys sample rate
    

   
    
    
    %% outputs
    
    
    %figs
    plot(f.fluo(1, :), 'Color',[0.9 0 0])
    box off
    %title('Background subtracted fluorescence', 'FontSize',14)
    ax=gca;
    ax.FontSize = 14;
    xlabel(figs.xtitle, 'FontSize',18)
    ylabel('fluorescence', 'FontSize',18)
    savefig(fullfile(path.figs,strcat(path.name,'_fluorescence.fig')))
    
    
 
    %saving workspace%
    save(fullfile(path.vars,strcat(path.name,'_workspace.mat')));
    
    close all


    
    %% repeate? %%
    choose = questdlg('process annother recording trial with these parameters?  ',...
        'continue or complete','yes','no','no');
    switch choose
        case 'yes'
            disp(' ');
            disp(' ');
            disp(' ');
            disp('select next recording directory');
            disp(' ');
            disp(' ');
            disp(' ');
            choose=-1;
        case 'no'
            disp('   analysis complete');
            choose=100;
    end
    
end

clear

