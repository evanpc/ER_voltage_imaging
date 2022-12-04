%% general notes %%




%FUNCTION:

%this script imports a single color multiplane tiff image, extracts the average fluorescence from windows of a designated width (xaxis) along a hand drawn ROI at intervals of a designated distance. Average background is subtracted from corresponding windows in a hand drawn background ROI. 

%this script cannot be done iteratively

%this script creates a variables and a figures directory in the home directory of the .tif image being analyzed. In these directories the workspace and fluorescence vs time image are saved, respectively. 

%the output variable for background subtracted fluorescence from each window is in theœrois.windows.f file.




%DATA REQUIREMENTS:

%a .tiff image is all

%FUNCTION REQUIREMENTS:

%this script also requires the readinim.m function




%INSTRUCTIONS:

%in general, follow instructions on UI boxes and check command line feedback. More explanation below:

%select directory where .tif image is stored
%change input parameters on GUI as needed
%draw a line along image that denotes the center of travel for the ROI. If needed you will be prompted to draw further distance.
%draw cell ROI completely encompassing center line.
%draw background ROI same width (in x dimension) as cell ROI, but over a background region of the image







%% Choosing images %%

%path to voltage image%
disp(' ');
disp(' ');
disp('load voltage image')
[ref.paths.vim.filnam,ref.paths.vim.dirpath, ~] = uigetfile('.tif'); %vimage
ref.paths.vim.impath=fullfile(ref.paths.vim.dirpath,ref.paths.vim.filnam);

%command line feedback
disp(' ');
disp(' ');
disp(strcat('voltage image is:    "',ref.paths.vim.filnam,'"'));



%% making directories and save paths %%
mkdir(fullfile(ref.paths.vim.dirpath,'figures'));
ref.paths.figs=fullfile(ref.paths.vim.dirpath,'figures');
mkdir(fullfile(ref.paths.vim.dirpath,'variables'));
ref.paths.vars=fullfile(ref.paths.vim.dirpath,'variables');



%% reading in images %%
tif.vim=readinim(ref.paths.vim.impath);



%% Window parameters %%




%imput/check parameters
tools.figs.dlg.prompt={'pixel size (microns)','window width (microns)',...
    'total distance to slide window (microns)','step size (microns)',...
    'duration of each frame (ms)','timescale for x axis (in ms)',...
    'title for xaxis'};
tools.figs.dlg.title='modify analysis parameters as needed';
tools.figs.dlg.dims=[1 70];
tools.figs.dlg.definput={'0.108', '5', '50', '5', '100', '1000', 'seconds'};%preset but changable parameters
tools.paras=inputdlg(tools.figs.dlg.prompt,tools.figs.dlg.title,...
    tools.figs.dlg.dims,tools.figs.dlg.definput);

%converting
tools.figs.xaxtit=char(tools.paras(7));
tools.paras=str2double(tools.paras);



%% hand drawn dendrite line for sliding rois %%


%tracing dendrite for roi info%
tools.choose=-1;
tools.figs.tit1=sprintf('begin where windows start');
tools.figs.tit2=sprintf('draw line along dendrite region to be analyzed');
figure('Color','white')

%quality control loop to ensure dend line is correct
while tools.choose<0
    imshow(tif.vim.ref);
    title(tools.figs.tit2, 'FontSize', 14);
    [rois.dendline.lnpixsxy(:,1),rois.dendline.lnpixsxy(:,2),~,...
        rois.dendline.lnsegsxy(:,1),rois.dendline.lnsegsxy(:,2)]=improfile;

    %visual check for dend line plot
    imshow(tif.vim.ref);
    hold on
    plot(rois.dendline.lnpixsxy(:,1),rois.dendline.lnpixsxy(:,2),...
        'Color','r','LineWidth',2);
    hold off


    %checks%
    tools.choose = questdlg('does dendrite line look correct?  ',...
        'check dend line!','yes','no','yes');
    switch tools.choose
        case 'yes'


            %making sure dend line is long enough
            if size(rois.dendline.lnpixsxy,1)<... dend line after max <
                    (tools.paras(3)+... total length/window size
                    (tools.paras(2)*2))... 1/2 a window size
                    /tools.paras(1) % converting microns to pix

                disp(' ');
                disp(' ');
                disp('   extend dend line for set bound');
                tools.figs.tit2=sprintf(' redraw dend line further');

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

    plot(rois.dendline.lnpixsxy(1:...
        (tools.paras(3)+(2*tools.paras(2)))/... total slide length plus two windows to be safe
        tools.paras(1),1),... converting to pix
        rois.dendline.lnpixsxy(1:...
        (tools.paras(3)+(2*tools.paras(2)))/... total slide length plus two windows to be safe
        tools.paras(1),2),'Color','r','LineWidth',2);

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

    plot(rois.dendline.lnpixsxy(1:...
        (tools.paras(3)+(2*tools.paras(2)))/... total slide length plus two windows to be safe
        tools.paras(1),1),... converting to pix
        rois.dendline.lnpixsxy(1:...
        (tools.paras(3)+(2*tools.paras(2)))/... total slide length plus two windows to be safe
        tools.paras(1),2),'Color','r','LineWidth',2);

    visboundaries(rois.totdend.masks,'Color',[0.6 0 0],'LineWidth',2);

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
    'pixdistst',zeros(floor(tools.paras(3)/...
    tools.paras(4)),1),...
    'prest',zeros(floor(tools.paras(3)/...
    tools.paras(4)),1),...
    'postst',zeros(floor(tools.paras(3)/...
    tools.paras(4)),1),...
    'thetast',zeros(floor(tools.paras(3)/...
    tools.paras(4)),1),...
    'xst',zeros(floor(tools.paras(3)/...
    tools.paras(4)),1),...
    'pixdistend',zeros(floor(tools.paras(3)/...
    tools.paras(4)),1),...
    'preend',zeros(floor(tools.paras(3)/...
    tools.paras(4)),1),...
    'postend',zeros(floor(tools.paras(3)/...
    tools.paras(4)),1),...
    'thetaend',zeros(floor(tools.paras(3)/...
    tools.paras(4)),1),...
    'xend',zeros(floor(tools.paras(3)/...
    tools.paras(4)),1));


%populating it is direction dependent
if rois.dendline.lnpixsxy(end,1)<rois.dendline.lnpixsxy(1,1)

    for i=1:floor(tools.paras(3)/tools.paras(4))


        %finding start bound x%
        rois.dendline.winpts.pixdistst(i)=...
            (tools.paras(2)/tools.paras(1))-...
            ((0.5*tools.paras(2))/tools.paras(1))+...
            ((tools.paras(4)/tools.paras(1))*(i-1));

        rois.dendline.winpts.prest(i)=find(rois.dendline.lnsegsxy(:,4)<=...
            rois.dendline.winpts.pixdistst(i),1,'last');

        rois.dendline.winpts.postst(i)=find(rois.dendline.lnsegsxy(:,4)>...
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
        rois.dendline.winpts.pixdistend(i)=...
            (tools.paras(2)/tools.paras(1))+...
            ((0.5*tools.paras(2))/tools.paras(1))+...
            ((tools.paras(4)/tools.paras(1))*(i-1));

        rois.dendline.winpts.preend(i)=find(rois.dendline.lnsegsxy(:,4)<=...
            rois.dendline.winpts.pixdistend(i),1,'last');

        rois.dendline.winpts.postend(i)=find(rois.dendline.lnsegsxy(:,4)>...
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


    for i=1:floor(tools.paras(3)/tools.paras(4))


        %finding start bound x%
        rois.dendline.winpts.pixdistst(i)=...
            (tools.paras(2)/tools.paras(1))-...
            ((0.5*tools.paras(2))/tools.paras(1))+...
            ((tools.paras(4)/tools.paras(1))*(i-1));

        rois.dendline.winpts.prest(i)=find(rois.dendline.lnsegsxy(:,4)<=...
            rois.dendline.winpts.pixdistst(i),1,'last');

        rois.dendline.winpts.postst(i)=find(rois.dendline.lnsegsxy(:,4)>...
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
        rois.dendline.winpts.pixdistend(i)=...
            (tools.paras(2)/tools.paras(1))+...
            ((0.5*tools.paras(2))/tools.paras(1))+...
            ((tools.paras(4)/tools.paras(1))*(i-1));

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
rois.windows(floor(tools.paras(3)/tools.paras(4)))=struct(...
    'mask',[],'fluo',[],'f',[],'compiled',[],'curvestats',[]);
rois.bkgrnd(size(rois.windows,2))=struct('mask',[],'fluo',[]);

%making restricted rois (direction dependent)%
if rois.dendline.lnpixsxy(end,1)<rois.dendline.lnpixsxy(1,1)
    for i=1:floor(tools.paras(3)/tools.paras(4))


        %dend windows first
        rois.windows(i).mask=rois.totdend.masks;
        rois.windows(i).mask(:,1:rois.dendline.winpts.xend(i))=0;
        rois.windows(i).mask(:,rois.dendline.winpts.xst(i):end)=0;

        rois.windows(i).fluo=zeros(1,tif.vim.xyz(3));
        rois.windows(i).f=zeros(1,tif.vim.xyz(3));

        %then bkgrnd
        rois.bkgrnd(i).mask=rois.totdend.bkmasks;
        rois.bkgrnd(i).mask(:,1:rois.dendline.winpts.xend(i))=0;
        rois.bkgrnd(i).mask(:,rois.dendline.winpts.xst(i):end)=0;

        rois.bkgrnd(i).fluo=zeros(1,tif.vim.xyz(3));
    end
else
    for i=1:floor(tools.paras(3)/tools.paras(4))


        %dend windows first
        rois.windows(i).mask=rois.totdend.masks;
        rois.windows(i).mask(:,rois.dendline.winpts.xend(i):end)=0;
        rois.windows(i).mask(:,1:rois.dendline.winpts.xst(i))=0;

        rois.windows(i).fluo=zeros(1,tif.vim.xyz(3));
        rois.windows(i).f=zeros(1,tif.vim.xyz(3));


        %then bkgrnd
        rois.bkgrnd(i).mask=rois.totdend.bkmasks;
        rois.bkgrnd(i).mask(:,rois.dendline.winpts.xend(i):end)=0;
        rois.bkgrnd(i).mask(:,1:rois.dendline.winpts.xst(i))=0;

        rois.bkgrnd(i).fluo=zeros(1,tif.vim.xyz(3));
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
        rois.windows(j).f=rois.windows(j).fluo-rois.bkgrnd(j).fluo; %net fluorescence
    end
end





%% figures %%




%full recording all overlapping%
tools.figs.fullxax=linspace((1*tools.paras(5)/...
    tools.paras(6)),(tif.vim.xyz(3)*tools.paras(5)/...
    tools.paras(6)),tif.vim.xyz(3));

plot(tools.figs.fullxax,rois.windows(1).f)
hold on

for j=2:size(rois.windows,2)
    plot(tools.figs.fullxax,rois.windows(j).f)
end

box off
xlabel(tools.figs.xaxtit, 'FontSize',18)
ylabel('dFoF', 'FontSize',18)
savefig(fullfile(ref.paths.figs,'full_fluo_all.fig'));


%total dendritic roi%
imshow(tif.vim.ref);
hold on
visboundaries(rois.totdend.masks,...
    'Color','r','LineWidth',2);
visboundaries(rois.totdend.bkmasks,...
    'Color','c','LineWidth',2);
hold off
savefig(fullfile(ref.paths.figs,'roi_totdend.fig'));



%close all



%save(fullfile(ref.paths.vars,strcat('slide_win_workspace_',extractBefore(ref.paths.vim.filnam,'_MMStack.ome.tif'),'.mat')));



%clear




