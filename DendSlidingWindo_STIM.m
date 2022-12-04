%% general notes %%




%FUNCTION:

%this script imports two single color multiplane tiff images - one of 
% ASAP3ER in dendrites and one of Alexa-549 fluorescence indicating caffeine 
% plumes. This script extracts the average fluorescence from windows of a 
% designated width (set these parameters by hand in first section belwo) 
% along a hand drawn ROI at intervals of a designated distance from the 
% center of caffeine application (which is determined automatically from 
% the Alexa-594 image). Average background is subtracted from corresponding 
% windows in a hand drawn background ROI. 

%this script cannot be done iteratively

%this script creates a variables and figures directory in the home 
% directory of the .tif image of ASAP3ER in dendrites. In these directories 
% the workspace and fluorescence vs time image are saved, respectively. 

%the output variable for background subtracted fluorescence from each 
% window is in the rois.windows.f(n) file.




%DATA INPUT REQUIREMENTS:

%two .tiff images

%FUNCTION REQUIREMENTS:

%this script also requires the readinim.m, nrm_cmple.m, bksub_cmple.m, and curve_stat_set.m functions


%INSTRUCTIONS:

%in general, follow instructions on UI boxes and check command line 
% feedback. More explanation below:

%change parameters of experiment in sections below by hand (GUI interface to come)%
%select directory where .tif image is stored
%draw a line along image that denotes the center of travel for the ROI. If needed you will be prompted to draw further distance.%
%draw cell ROI completely encompassing center line.
%draw background ROI same width (in x dimension) as cell ROI, but over a background region of the image%


%%


%% ADJUST THESE PARAMETERS BY HAND FOR EACH IMAGE IF NEEDED %%

%%


%% Window parameters %%

%pixel size in microns
tools.parameters.windows(1)=0.108;
%window size in microns
tools.parameters.windows(2)=6;
%total distance to slide window in microns
tools.parameters.windows(3)=50;
%step size in microns
tools.parameters.windows(4)=0.5;



%% Image parameters %%

%image duration in ms
tools.parameters.im(1)=100;
%scale for xaxis in ms
tools.parameters.im(2)=1000;
%corrisponding time unit for xaxis
tools.figs.xaxtit='seconds';
%expected number of images for double checking
tools.parameters.im(3)=1190;



%% Stim parameters %%

%time of stim pulses in ms
tools.parameters.stimtimes=[15000,45000,75000];
%time viewed before in ms
tools.parameters.stim(1)=4900;
%time viewed after in ms
tools.parameters.stim(2)=25000;
%deflection (negative is -1, pos is 1)
tools.parameters.stim(3)=-1;
%reasonable peak window in ms
tools.parameters.stim(4)=1000;


%%

warning('off', 'all')



%% Choosing images %%

%path to voltage image%
disp(' ');
disp(' ');
disp('load voltage image')
[ref.paths.vim.filnam,ref.paths.vim.dirpath, ~] = uigetfile('.tif'); %vimage
ref.paths.vim.impath=fullfile(ref.paths.vim.dirpath,ref.paths.vim.filnam);

%path to 594 reference image
disp(' ');
disp(' ');
disp('load 594 reference image')
[ref.paths.im594.filnam,ref.paths.im594.dirpath, ~]=...
    uigetfile(fullfile(ref.paths.vim.dirpath,'*.tif')); %594 image
ref.paths.im594.impath=fullfile(ref.paths.im594.dirpath,...
    ref.paths.im594.filnam);

%command line feedback
disp(' ');
disp(' ');
disp(strcat('voltage image is:    "',ref.paths.vim.filnam,'"'));
disp(strcat('594 reference image is:    "',ref.paths.im594.filnam,'"'));



%% making directories and save paths %%
mkdir(fullfile(ref.paths.vim.dirpath,'figures'));
ref.paths.figs=fullfile(ref.paths.vim.dirpath,'figures');
mkdir(fullfile(ref.paths.vim.dirpath,'variables'));
ref.paths.vars=fullfile(ref.paths.vim.dirpath,'variables');



%% reading in images %%
tif.vim=readinim(ref.paths.vim.impath);
tif.im594=readinim(ref.paths.im594.impath);



%% hand drawn dendrite line for sliding rois %%


%tracing dendrite for roi info%
tools.choose=-1;
tools.figs.tit1=sprintf('reference 594 application heatmap');
tools.figs.tit2=sprintf('draw line along dendrite region to be analyzed');
figure('Color','white')

%quality control loop to ensure dend line is correct
while tools.choose<0
    subplot(2,1,1);
    imshow(tif.im594.ref);
    colormap('jet') %change jet if you want to change heatmap type
    title(tools.figs.tit1, 'FontSize', 14);
    subplot(2,1,2);
    imshow(tif.vim.ref);
    title(tools.figs.tit2, 'FontSize', 14);
    [rois.dendline.lnpixsxy(:,1),rois.dendline.lnpixsxy(:,2),~,...
        rois.dendline.lnsegsxy(:,1),rois.dendline.lnsegsxy(:,2)]=improfile;
    
    %visual check for dend line plot
    subplot(2,1,2);
    imshow(tif.vim.ref);
    hold on
    plot(rois.dendline.lnpixsxy(:,1),rois.dendline.lnpixsxy(:,2),...
        'Color','r','LineWidth',2);
    
    
    %checks%
    tools.choose = questdlg('does dendrite line look correct?  ',...
        'check dend line!','yes','no','yes');
    switch tools.choose
        case 'yes'
            
            %finding start point from im594 ref centroid
            rois.dendline.prof594=improfile(tif.im594.ref,...
                rois.dendline.lnsegsxy(:,1),rois.dendline.lnsegsxy(:,2));
            rois.dendline.lnmax.cell=round(mean(find(rois.dendline.prof594==...
                max(rois.dendline.prof594))));
            
            %making sure dend line is long enough
            if size(rois.dendline.prof594)-rois.dendline.lnmax.cell<... dend line after max <
                    (tools.parameters.windows(3)+... total length/window size
                    (tools.parameters.windows(2)/2))... 1/2 a window size
                    /tools.parameters.windows(1) % converting microns to pix
                if rois.dendline.lnmax.cell<... start point
                        (tools.parameters.windows(2)/2)... half window size
                        /tools.parameters.windows(1) % converting microns to pix
                    disp(' ');
                    disp(' ');
                    disp('   start dend line earlier and extend longer');
                    tools.figs.tit1=sprintf(' start dend line earlier');
                else
                    disp(' ');
                    disp(' ');
                    disp('   extend dend line further for set bound');
                    tools.figs.tit2=sprintf(' redraw dend line further');
                end
                tools.choose=-1;
                rois.dendline=rmfield(rois.dendline,'lnpixsxy');
                rois.dendline=rmfield(rois.dendline,'lnsegsxy');
            elseif rois.dendline.lnmax.cell<... start point
                    (tools.parameters.windows(2)/2)... half window size
                    /tools.parameters.windows(1) % converting microns to pix
                disp(' ');
                disp(' ');
                disp('   start dend line earlier');
                tools.figs.tit1=sprintf(' start dend line earlier');
                tools.choose=-1;
                rois.dendline=rmfield(rois.dendline,'lnpixsxy');
                rois.dendline=rmfield(rois.dendline,'lnsegsxy');
            else, tools.choose=1;
                disp(' ');
                disp(' ');
                disp('   dend line set')
            end
        case 'no'
            tools.figs.tit2=sprintf(' redraw dend line');
            tools.choose =-1;
            rois.dendline=rmfield(rois.dendline,'lnpixsxy');
            rois.dendline=rmfield(rois.dendline,'lnsegsxy');
    end
end
savefig(fullfile(ref.paths.figs,'roi_dendline.fig'));
close all



%% drawing roi frame for entire dendrite region %%


%drawing master roi for windows to be built from%
tools.figs.tit1=sprintf('draw roi around dendrite (at least as far as line)');
figure('Color','white')

while tools.choose<2
    
    %creating image with line ref for minimal roi length
    imshow(tif.vim.ref);
    title(tools.figs.tit1, 'FontSize', 14);
    hold on
    plot(rois.dendline.lnpixsxy((rois.dendline.lnmax.cell-...
        (tools.parameters.windows(2)/2)/tools.parameters.windows(1)):...
        ((tools.parameters.windows(3)+(tools.parameters.windows(2)/2))/...
        tools.parameters.windows(1))+rois.dendline.lnmax.cell,1),...
        rois.dendline.lnpixsxy((rois.dendline.lnmax.cell-...
        (tools.parameters.windows(2)/2)/tools.parameters.windows(1)):...
        ((tools.parameters.windows(3)+(tools.parameters.windows(2)/2))/...
        tools.parameters.windows(1))+rois.dendline.lnmax.cell,2),...
        'Color','r','LineWidth',2);
    hold off
    rois.totdend.frm = imfreehand;
    setColor(rois.totdend.frm,'r');
    
    %NOTE: points sampled by improfile (linpixsxy) corrdinates are used to
    %generate line shown on plot, which is a slight overestimate but is
    %fine for this approximation.  More accuracy is used later for
    %generating window bounds.
    
    tools.choose = questdlg('cell roi correct?  ',...
        'check cell roi!','yes','no','yes');
    switch tools.choose
        case 'yes'
            disp(' ');
            disp(' ');
            disp('   full dendrite roi set')
            tools.choose = 2;
        case 'no'
            tools.figs.tit2=sprintf(' redraw dendrite roi');
            tools.choose =-1;
    end
end
rois.totdend.masks=createMask(rois.totdend.frm);
savefig(fullfile(ref.paths.figs,'roi_fulldend.fig'));
close all

%making backgound mask
tools.figs.tit1=sprintf('draw roi for background (at least as wide as line)');
figure('Color','white')
while tools.choose<3
    imshow(tif.vim.ref);
    title(tools.figs.tit1, 'FontSize', 14);
    hold on
    plot(rois.dendline.lnpixsxy((rois.dendline.lnmax.cell-...
        (tools.parameters.windows(2)/2)/tools.parameters.windows(1)):...
        ((tools.parameters.windows(3)+(tools.parameters.windows(2)/2))/...
        tools.parameters.windows(1))+rois.dendline.lnmax.cell,1),...
        rois.dendline.lnpixsxy((rois.dendline.lnmax.cell-...
        (tools.parameters.windows(2)/2)/tools.parameters.windows(1)):...
        ((tools.parameters.windows(3)+(tools.parameters.windows(2)/2))/...
        tools.parameters.windows(1))+rois.dendline.lnmax.cell,2),...
        'Color','r','LineWidth',2);
    hold off
    rois.totdend.bkfrm = imfreehand;
    setColor(rois.totdend.bkfrm,'c');
    
    tools.choose = questdlg('background roi correct?  ',...
        'check background roi!','yes','no','yes');
    switch tools.choose
        case 'yes'
            tools.choose = 3;
            disp(' ');
            disp(' ');
            disp('   background set')
        case 'no'
            tools.figs.tit1=('redraw background roi');
            tools.choose =-1;
    end
end
rois.totdend.bkmasks = createMask(rois.totdend.bkfrm);%this var may need to not exist
savefig(fullfile(ref.paths.figs,'roi_bkgrnd.fig'));



%% Calculating window start and end bounds from dendrite line distances%%

%pythag for pix distances between points and along line for every point
rois.dendline.lnsegsxy(:,3:4)=0;
for i=2:size(rois.dendline.lnsegsxy,1)
    rois.dendline.lnsegsxy(i,3)=sqrt((abs(rois.dendline.lnsegsxy(i,1)-...
        rois.dendline.lnsegsxy(i-1,1)).^2)+(abs(rois.dendline.lnsegsxy(i,2)-...
        rois.dendline.lnsegsxy(i-1,2)).^2));
    rois.dendline.lnsegsxy(i,4)=rois.dendline.lnsegsxy(i,3)+...
        rois.dendline.lnsegsxy(i-1,4);
end


%generating structure with accurate x bounds for roi masks%
rois.dendline.winpts=struct(...
    'pixdistst',zeros(floor(tools.parameters.windows(3)/...
    tools.parameters.windows(4)),1),...
    'prest',zeros(floor(tools.parameters.windows(3)/...
    tools.parameters.windows(4)),1),...
    'postst',zeros(floor(tools.parameters.windows(3)/...
    tools.parameters.windows(4)),1),...
    'thetast',zeros(floor(tools.parameters.windows(3)/...
    tools.parameters.windows(4)),1),...
    'xst',zeros(floor(tools.parameters.windows(3)/...
    tools.parameters.windows(4)),1),...
    'pixdistend',zeros(floor(tools.parameters.windows(3)/...
    tools.parameters.windows(4)),1),...
    'preend',zeros(floor(tools.parameters.windows(3)/...
    tools.parameters.windows(4)),1),...
    'postend',zeros(floor(tools.parameters.windows(3)/...
    tools.parameters.windows(4)),1),...
    'thetaend',zeros(floor(tools.parameters.windows(3)/...
    tools.parameters.windows(4)),1),...
    'xend',zeros(floor(tools.parameters.windows(3)/...
    tools.parameters.windows(4)),1));


%populating it is direction dependent
if rois.dendline.lnpixsxy(end,1)<rois.dendline.lnpixsxy(1,1)
    
    
    %determining distance of application midpoint%
    rois.dendline.lnmax.preseg=find(rois.dendline.lnsegsxy(:,1)>...
        rois.dendline.lnpixsxy(rois.dendline.lnmax.cell,1),1,'last');
    rois.dendline.lnmax.totpixdist=...
        rois.dendline.lnsegsxy(rois.dendline.lnmax.preseg,4)+...
        (sqrt(((rois.dendline.lnsegsxy(rois.dendline.lnmax.preseg,1)-...
        rois.dendline.lnpixsxy(rois.dendline.lnmax.cell,1)).^2)+...
        ((rois.dendline.lnsegsxy(rois.dendline.lnmax.preseg,2)-...
        rois.dendline.lnpixsxy(rois.dendline.lnmax.cell,2)).^2)));
    
    for i=1:floor(tools.parameters.windows(3)/tools.parameters.windows(4))
        
        
        %finding start bound x%
        rois.dendline.winpts.pixdistst(i)=rois.dendline.lnmax.totpixdist-...
            ((0.5*tools.parameters.windows(2))/tools.parameters.windows(1))+...
            ((tools.parameters.windows(4)/tools.parameters.windows(1))*(i-1));
        rois.dendline.winpts.prest(i)=find(rois.dendline.lnsegsxy(:,4)<...
            rois.dendline.winpts.pixdistst(i),1,'last');
        rois.dendline.winpts.postst(i)=find(rois.dendline.lnsegsxy(:,4)>=...
            rois.dendline.winpts.pixdistst(i),1,'first');
        rois.dendline.winpts.thetast(i)=atan(...
            (rois.dendline.lnsegsxy(rois.dendline.winpts.postst(i),2)-...
            rois.dendline.lnsegsxy(rois.dendline.winpts.prest(i),2))/...
            (rois.dendline.lnsegsxy(rois.dendline.winpts.postst(i),1)-...
            rois.dendline.lnsegsxy(rois.dendline.winpts.prest(i),1)));
        rois.dendline.winpts.xst(i)=...
            rois.dendline.lnsegsxy(rois.dendline.winpts.prest(i),1)-... %change this - for oposite direction
            (round(cos(rois.dendline.winpts.thetast(i))*...
            (rois.dendline.winpts.pixdistst(i)-...
            rois.dendline.lnsegsxy(rois.dendline.winpts.prest(i),4))));
        
        
        %finding end bound x%
        rois.dendline.winpts.pixdistend(i)=rois.dendline.lnmax.totpixdist+...
            ((0.5*tools.parameters.windows(2))/tools.parameters.windows(1))+...
            ((tools.parameters.windows(4)/tools.parameters.windows(1))*(i-1));
        rois.dendline.winpts.preend(i)=find(rois.dendline.lnsegsxy(:,4)<...
            rois.dendline.winpts.pixdistend(i),1,'last');
        rois.dendline.winpts.postend(i)=find(rois.dendline.lnsegsxy(:,4)>=...
            rois.dendline.winpts.pixdistend(i),1,'first');
        rois.dendline.winpts.thetaend(i)=atan(...
            (rois.dendline.lnsegsxy(rois.dendline.winpts.postend(i),2)-...
            rois.dendline.lnsegsxy(rois.dendline.winpts.preend(i),2))/...
            (rois.dendline.lnsegsxy(rois.dendline.winpts.postend(i),1)-...
            rois.dendline.lnsegsxy(rois.dendline.winpts.preend(i),1)));
        rois.dendline.winpts.xend(i)=...
            rois.dendline.lnsegsxy(rois.dendline.winpts.preend(i),1)-... %change this - for oposite direction
            (round(cos(rois.dendline.winpts.thetaend(i))*...
            (rois.dendline.winpts.pixdistend(i)-...
            rois.dendline.lnsegsxy(rois.dendline.winpts.preend(i),4))));
    end
    
else
    
    
    %determining distance of application midpoint%
    rois.dendline.lnmax.preseg=find(rois.dendline.lnsegsxy(:,1)<...
        rois.dendline.lnpixsxy(rois.dendline.lnmax.cell,1),1,'last');
    rois.dendline.lnmax.totpixdist=...
        rois.dendline.lnsegsxy(rois.dendline.lnmax.preseg,4)+...
        (sqrt(((rois.dendline.lnpixsxy(rois.dendline.lnmax.cell,1)-...
        rois.dendline.lnsegsxy(rois.dendline.lnmax.preseg,1)).^2)+...
        ((rois.dendline.lnpixsxy(rois.dendline.lnmax.cell,2)-...
        rois.dendline.lnsegsxy(rois.dendline.lnmax.preseg,2)).^2)));
    
    for i=1:floor(tools.parameters.windows(3)/tools.parameters.windows(4))
        
        
        %finding start bound x%
        rois.dendline.winpts.pixdistst(i)=rois.dendline.lnmax.totpixdist-...
            ((0.5*tools.parameters.windows(2))/tools.parameters.windows(1))+...
            ((tools.parameters.windows(4)/tools.parameters.windows(1))*(i-1));
        rois.dendline.winpts.prest(i)=find(rois.dendline.lnsegsxy(:,4)<...
            rois.dendline.winpts.pixdistst(i),1,'last');
        rois.dendline.winpts.postst(i)=find(rois.dendline.lnsegsxy(:,4)>=...
            rois.dendline.winpts.pixdistst(i),1,'first');
        rois.dendline.winpts.thetast(i)=atan(...
            (rois.dendline.lnsegsxy(rois.dendline.winpts.postst(i),2)-...
            rois.dendline.lnsegsxy(rois.dendline.winpts.prest(i),2))/...
            (rois.dendline.lnsegsxy(rois.dendline.winpts.postst(i),1)-...
            rois.dendline.lnsegsxy(rois.dendline.winpts.prest(i),1)));
        rois.dendline.winpts.xst(i)=round(...
            rois.dendline.lnsegsxy(rois.dendline.winpts.prest(i),1)+... %change this + for oposite direction
            (round(cos(rois.dendline.winpts.thetast(i))*...
            (rois.dendline.winpts.pixdistst(i)-...
            rois.dendline.lnsegsxy(rois.dendline.winpts.prest(i),4)))));
        
        
        %finding end bound x%
        rois.dendline.winpts.pixdistend(i)=rois.dendline.lnmax.totpixdist+...
            ((0.5*tools.parameters.windows(2))/tools.parameters.windows(1))+...
            ((tools.parameters.windows(4)/tools.parameters.windows(1))*(i-1));
        rois.dendline.winpts.preend(i)=find(rois.dendline.lnsegsxy(:,4)<...
            rois.dendline.winpts.pixdistend(i),1,'last');
        rois.dendline.winpts.postend(i)=find(rois.dendline.lnsegsxy(:,4)>=...
            rois.dendline.winpts.pixdistend(i),1,'first');
        rois.dendline.winpts.thetaend(i)=atan(...
            (rois.dendline.lnsegsxy(rois.dendline.winpts.postend(i),2)-...
            rois.dendline.lnsegsxy(rois.dendline.winpts.preend(i),2))/...
            (rois.dendline.lnsegsxy(rois.dendline.winpts.postend(i),1)-...
            rois.dendline.lnsegsxy(rois.dendline.winpts.preend(i),1)));
        rois.dendline.winpts.xend(i)=round(...
            rois.dendline.lnsegsxy(rois.dendline.winpts.preend(i),1)+... %change this + for oposite direction
            (round(cos(rois.dendline.winpts.thetaend(i))*...
            (rois.dendline.winpts.pixdistend(i)-...
            rois.dendline.lnsegsxy(rois.dendline.winpts.preend(i),4)))));
    end
end



%% making windows %%


%generating structure with accurate x bounds for roi masks%
rois.windows(floor(tools.parameters.windows(3)/tools.parameters.windows(4)))=struct(...
    'mask',[],'fluo',[],'f',[],'compiled',[],'curvestats',[],...
    'raw594',[],'compiled594',[],'curvestats594',[]);
rois.bkgrnd(size(rois.windows,2))=struct('mask',[],'fluo',[]);

%making restricted rois (direction dependent)%
if rois.dendline.lnpixsxy(end,1)<rois.dendline.lnpixsxy(1,1)
    for i=1:floor(tools.parameters.windows(3)/tools.parameters.windows(4))
        
        
        %dend windows first
        rois.windows(i).mask=rois.totdend.masks;
        rois.windows(i).mask(:,1:rois.dendline.winpts.xend(i))=0;
        rois.windows(i).mask(:,rois.dendline.winpts.xst(i):end)=0;
        
        rois.windows(i).fluo=zeros(1,tif.vim.xyz(3));
        rois.windows(i).f=zeros(1,tif.vim.xyz(3));
        rois.windows(i).raw594=zeros(1,tif.im594.xyz(3));
        
        %then bkgrnd
                rois.bkgrnd(i).mask=rois.totdend.bkmasks;
        rois.bkgrnd(i).mask(:,1:rois.dendline.winpts.xend(i))=0;
        rois.bkgrnd(i).mask(:,rois.dendline.winpts.xst(i):end)=0;
        
        rois.bkgrnd(i).fluo=zeros(1,tif.im594.xyz(3));
    end
else
    for i=1:floor(tools.parameters.windows(3)/tools.parameters.windows(4))
        
        
        %dend windows first
        rois.windows(i).mask=rois.totdend.masks;
        rois.windows(i).mask(:,rois.dendline.winpts.xend(i):end)=0;
        rois.windows(i).mask(:,1:rois.dendline.winpts.xst(i))=0;
        
        rois.windows(i).fluo=zeros(1,tif.vim.xyz(3));
        rois.windows(i).f=zeros(1,tif.vim.xyz(3));
        rois.windows(i).raw594=zeros(1,tif.im594.xyz(3));
        
        
        %then bkgrnd
                rois.bkgrnd(i).mask=rois.totdend.bkmasks;
        rois.bkgrnd(i).mask(:,rois.dendline.winpts.xend(i):end)=0;
        rois.bkgrnd(i).mask(:,1:rois.dendline.winpts.xst(i))=0;
        
        rois.bkgrnd(i).fluo=zeros(1,tif.im594.xyz(3));
    end
end



%% reading fluorescence values %%


%reading average fluo from respective rois for each frame%
for i=1:size(tif.vim.reads,3)
    
    %pulling in one z plane
    tools.temp=tif.vim.reads(:,:,i);
    
    %reading rois from each z plane
    for j=1:size(rois.windows,2)
        rois.windows(j).fluo(i)=mean(tools.temp(rois.windows(j).mask)); %mean of roi on v image
            rois.bkgrnd(j).fluo(i)=mean(tools.temp(rois.bkgrnd(j).mask)); %mean of roi on image
    end
end


%reading average 594 raw fluo at each roi for each frame
for i=1:size(tif.im594.reads,3)
    
    %pulling in one z plane
    tools.temp=tif.im594.reads(:,:,i);
    
    %reading rois from each z plane
    for j=1:size(rois.windows,2)
        rois.windows(j).raw594(i)=mean(tools.temp(rois.windows(j).mask));  %mean of roi on 594 image
    end
    
end


%double checking size and making them identical if needed
if size(tif.im594.reads,3)>size(tif.vim.reads,3)  
    for j=1:size(rois.windows,2)
        temp=(ones(1,(size(tif.im594.reads,3)-size(tif.vim.reads,3)))).*...
            rois.windows(j).fluo(1);
        rois.windows(j).fluo=[temp,rois.windows(j).fluo];
    end
    disp(' ');
    disp(' ');
    disp('voltage image is short');   
elseif size(tif.vim.reads,3)>size(tif.im594.reads,3)  
    for j=1:size(rois.windows,2)
        temp=(ones(1,(size(tif.vim.reads,3)-size(tif.im594.reads,3)))).*...
            rois.windows(j).raw594(1);
        rois.windows(j).raw594=[temp,rois.windows(j).raw594];
    end
    disp(' ');
    disp(' ');
    disp('594 image is short'); 
else
    disp(' ');
    disp(' ');
    disp('voltage image and 594 image are equal lengths'); 
end



%% compiling into averages and stats %%


%making summary variable as a function of distance
fxnodist.xax=(linspace(0,(size(rois.windows,2)-1),size(rois.windows,2))*...
    tools.parameters.windows(4));
fxnodist.amp=zeros(1,size(rois.windows,2));
fxnodist.tmax=zeros(1,size(rois.windows,2));
fxnodist.amp594=zeros(1,size(rois.windows,2));
fxnodist.tmax594=zeros(1,size(rois.windows,2));


%subtracting background, compiling resonses, generating stats%
for j=1:size(rois.windows,2)
    
    %backgrounded fluorescence
    rois.windows(j).f=rois.windows(j).fluo-rois.bkgrnd(j).fluo;
    
    %compiling average response curve from stims
    rois.windows(j).compiled=nrm_cmple(rois.windows(j).f,...
        (tools.parameters.stimtimes/tools.parameters.im(1)),...
        (tools.parameters.stim(1)/tools.parameters.im(1)),...
        (tools.parameters.stim(2)/tools.parameters.im(1)));
    
    rois.windows(j).compiled594=bksub_cmple(rois.windows(j).raw594,...
        (tools.parameters.stimtimes/tools.parameters.im(1)),...
        (tools.parameters.stim(1)/tools.parameters.im(1)),...
        (tools.parameters.stim(2)/tools.parameters.im(1)));
    
    %calculating curve stats
    rois.windows(j).curvestats=curve_stat_set(...
        (tools.parameters.stim(3)*rois.windows(j).compiled.stats.ave),...
        (tools.parameters.stim(1)/tools.parameters.im(1)),...
        (tools.parameters.stim(4)/tools.parameters.im(1)));
    
    rois.windows(j).curvestats594=curve_stat_set(...
        rois.windows(j).compiled594.stats.ave,...
        (tools.parameters.stim(1)/tools.parameters.im(1)),...
        (tools.parameters.stim(4)/tools.parameters.im(1)));
    
    %adding to summary variable
    fxnodist.amp(j)=rois.windows(j).curvestats.max(1);
    fxnodist.tmax(j)=rois.windows(j).curvestats.max(2)-(1+...
        tools.parameters.stim(1)/tools.parameters.im(1));
    fxnodist.amp594(j)=rois.windows(j).curvestats594.max(1);
    fxnodist.tmax594(j)=rois.windows(j).curvestats594.max(2)-(1+...
        tools.parameters.stim(1)/tools.parameters.im(1));
    
end

%estimating lambda from points plotted
[fxnodist.maxamp(1),fxnodist.maxamp(2)]=max(...
    fxnodist.amp);



%% figures %%


%full recording figures%
tools.figs.fullxax=linspace((tools.parameters.im(1)/...
    tools.parameters.im(2)),((size(rois.windows(1).f,2))...
    *tools.parameters.im(1)/tools.parameters.im(2)),...
    size(rois.windows(1).f,2));

plot(tools.figs.fullxax,rois.windows(fxnodist.maxamp(2)).f,...
    'Color',[0 0.6 0])
box off
title('Max fluorescence response', 'FontSize',18)
xlabel(tools.figs.xaxtit, 'FontSize',18)
ylabel('dFoF', 'FontSize',18)
savefig(fullfile(ref.paths.figs,'full_fluo_maxamp.fig'));


%rois
imshow(tif.vim.ref);
hold on
visboundaries(rois.totdend.masks,...
    'Color','r','LineWidth',2);
visboundaries(rois.totdend.bkmasks,...
    'Color','c','LineWidth',2);
hold off
title('max resoponse roi', 'FontSize', 18);
savefig(fullfile(ref.paths.figs,'roi_totdend.fig'));

imshow(tif.vim.ref);
hold on
visboundaries(rois.windows(fxnodist.maxamp(2)).mask,...
    'Color','r','LineWidth',2);
visboundaries(rois.bkgrnd(fxnodist.maxamp(2)).mask,...
    'Color','c','LineWidth',2);
hold off
title('max resoponse roi', 'FontSize', 18);
savefig(fullfile(ref.paths.figs,'roi_maxamp.fig'));



%paired images
tools.figs.left_color = [0 0.5 0];
tools.figs.right_color = [0.6 0.4 0.4];
tools.figs.bk_color = [1 1 1];


set(figure,'defaultAxesColorOrder',...
    [tools.figs.left_color; tools.figs.right_color],...
    'Color',tools.figs.bk_color);
yyaxis left
plot(fxnodist.xax,fxnodist.amp,'Color',[0 0.6 0],'LineWidth',1,...
    'Marker','o')
%axis([0 100000 -5 15])
title('Max amp and 594 fluorescence as a function of distance along dend', 'FontSize',18)
ylabel('dFoF max amplitude', 'FontSize',18)
hold on
yyaxis right
plot(fxnodist.xax,fxnodist.amp594,'Color',[0.55 0 0],'LineWidth',1,...
    'Marker','o')
box off
ylabel('594au', 'FontSize',18)
xlabel('distance to window center (um)', 'FontSize',18)
%axis([0 100000 -5 15])
hold off
savefig(fullfile(ref.paths.figs,'dist_vs_amp_both.fig'));


%set(figure,'defaultAxesColorOrder',...
    %[tools.figs.left_color; tools.figs.right_color],...
    %'Color',tools.figs.bk_color);
%yyaxis left
%plot(fxnodist.xax,fxnodist.tmax,'Color',[0 0.6 0],'LineWidth',1,...
    %'Marker','o')
%axis([0 100000 -5 15])
%ylabel('frames# after stim to max ASAP3er amp', 'FontSize',18)
%hold on
%yyaxis right
%plot(fxnodist.xax,fxnodist.tmax594,'Color',[0.55 0 0],'LineWidth',1,...
    %'Marker','o')
%box off
%ylabel('frame# after stim to max 594au', 'FontSize',18)
%xlabel('distance to window center (um)', 'FontSize',18)
%axis([0 100000 -5 15])
%hold off
%savefig(fullfile(ref.paths.figs,'dist_vs_tmax_both.fig'));


save(fullfile(ref.paths.vars,strcat('slide_win_workspace_',extractBefore(ref.paths.vim.filnam,'_MMStack.ome.tif'),'.mat')));

%tstim may need to have one subtracted. Currently this is done in the
%curve_stat_set fxn.


close all



%% Notes %%







