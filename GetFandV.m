%% general notes %%




%FUNCTION:

%this script imports a single color multiplane .tif image and a .mat 
%physiology variable. It extracts the average fluorescence from a hand 
%drawn ROI over each .tif plane (i.e. time), then subtracts the average 
%fluorescence from a hand drawn background ROI from the same plane. It then 
%linearly interpolates between each fluorescence point to fit the timescale 
%of physiology data, aligns fluorescence and voltage with respect to time, 
%and plots the two together. 

%this script can be done iteratively without reloading

%this script creates a variables and figures directory in the home 
%directory of the .tif image being analyzed. In these directories the 
%workspace and figures are saved, respectively. The output variables for 
%fluorescence and voltage aligned with each other are f.f and v.paired
%respectively. 




%REQUIREMENTS:

%a .tiff image and a .mat file that has physiology data in a field titled 
%“data”



%INSTRUCTIONS:

%in general, follow instructions on UI boxes and check command line feedback. More explanation below:

%select directory where .tif and .mat files are stored
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


%loading phys
v.name=dir(fullfile(path.expt,'*.mat'));
v.v=getfield(load(fullfile(v.name.folder,v.name.name)),...
    strrep(v.name.name,'.mat',''),'data');%making phys file name, loading phys file, loading it, accessing it, then the data, and saving data as v


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
disp(strcat('phys file is:    "',v.name.name,'"'));
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


figs.dlg.tit='modify experiment parameters as needed';
while choose<3
    figs.dlg.prompt={'phys sample period (ms)','image sample period (ms)',...
        strcat('ms before imaging (image frames=',num2str(tif.xyz(3)),...
        ', phys samples=',num2str(length(v.v)),')'),...
        strcat('ms after imaging (image frames=',num2str(tif.xyz(3)),...
        ', phys samples=',num2str(length(v.v)),')')};
    figs.dlg.title=figs.dlg.tit;
    figs.dlg.dims=[1 70];
    figs.dlg.definput={'0.1', '100', '500', '500'};%preset but changable parameters
    paras=str2double(inputdlg(figs.dlg.prompt,figs.dlg.title,...
        figs.dlg.dims,figs.dlg.definput));
    
    
    %these conversion varriables may not be necessary to create... but it condenses code below... maybe write this out in future%
    vn.pre=paras(3)/paras(1);
    vn.during=tif.xyz(3)*paras(2)/paras(1);
    vn.post=paras(4)/paras(1);
    vn.vpf=paras(2)/paras(1);
    
    
    %double checking that parameters line up
    if (vn.pre+vn.during+vn.post)==length(v.v)
        choose=3;
    else
        figs.dlg.tit='parameters are incorrect, re-enter';
        choose=-1;
    end
end



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
f.fluo=roi.fluo-bk.fluo;%this is fluo for experiment but not compatible withh phys sample rate


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


% matching f.fluo to phys time scale with linspace and adjusting x axis
f.f=[1:(vn.during-vn.vpf);zeros(1,vn.during-vn.vpf)];
f.f(1,:)=f.f(1,:).*tscale;   %changing row 1 acquisitions to time
for i=0:tif.xyz(3)-2 %adding appropriate point number in between frames
    f.f(2,(i*vn.vpf)+1:(i+1)*vn.vpf)=linspace(f.fluo(i+1),f.fluo(i+2),vn.vpf);
end


%including timescale to v and isolating only paired recording time
v.v=[1:length(v.v); v.v];
v.v(1,:)=v.v(1,:).*tscale;
v.paired=[f.f(1, :);v.v(2,(vn.pre+vn.vpf+1):(length(v.v)-vn.post))];



%% outputs




%figs
plot(f.f(1, :), f.f(2, :),'Color',[0.9 0 0])
box off
ax=gca;
ax.FontSize = 14;
xlabel(figs.xtitle, 'FontSize',18)
ylabel('fluorescence', 'FontSize',18)
savefig(fullfile(path.figs,strcat(path.name,'_fluorescence.fig')))

figs.left_color = [0 0.4470 0.7410];
figs.right_color = [0.85 0 0];
figs.bk_color = [1, 1, 1,];
set(figure,'defaultAxesColorOrder',[figs.left_color; figs.right_color],...
    'Color',figs.bk_color);
yyaxis right
plot(f.f(1, :), f.f(2, :),'Color',[0.85 0 0])
%axis([0 100000 -5 15])
title(path.name, 'FontSize',14)
ax1=gca;
ax1.FontSize = 14;
ylabel('Fluorescence', 'FontSize',12)
hold on
yyaxis left
plot(v.paired(1, :), v.paired(2, :))
box off
ax2=gca;
ax2.FontSize = 14;
ylabel('mV', 'FontSize',18)
xlabel(figs.xtitle, 'FontSize',18)
%axis([0 100000 -5 15])
hold off
savefig(fullfile(path.figs,strcat(path.name,'_paired.fig')))

close all


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
    
    
    %loading phys
    v.name=dir(fullfile(path.expt,'*.mat'));
    v.v=getfield(load(fullfile(v.name.folder,v.name.name)),...
        strrep(v.name.name,'.mat',''),'data');%making phys file name, loading phys file, loading it, accessing it, then the data, and saving data as v
    
    
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
    disp(strcat('phys file is:    "',v.name.name,'"'));
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
    
    
    vn.during=tif.xyz(3)*paras(2)/paras(1);
    
    if (vn.pre+vn.during+vn.post)==length(v.v)
        choose=choose+1;
    end
    
    figs.dlg.tit='new experiment parameters dont match! re-enter';
    figs.dlg.prompt={'phys sample period (ms)','image sample period (ms)',...
        strcat('ms before imaging (image frames=',num2str(tif.xyz(3)),...
        ', phys samples=',num2str(length(v.v)),')'),...
        strcat('ms after imaging (image frames=',num2str(tif.xyz(3)),...
        ', phys samples=',num2str(length(v.v)),')')};
    while choose<2
        paras = str2double(inputdlg(figs.dlg.prompt,figs.dlg.title,...
            figs.dlg.dims,figs.dlg.definput));
        
        %these conversion varriables may not be necessary to create... but it condenses code below... maybe write this out in future%
        vn.pre=paras(3)/paras(1);
        vn.during=tif.xyz(3)*paras(2)/paras(1);
        vn.post=paras(4)/paras(1);
        vn.vpf=paras(2)/paras(1);
        
        if (vn.pre+vn.during+vn.post)==length(v.v)
            choose=2;
        else
            figs.dlg.tit='parameters are incorrect, re-enter';
            choose=-1;
        end
    end
    
    
    
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
    
    
    % matching f.fluo to phys time scale with linspace and adjusting x axis
    f.f = [1:(vn.during-vn.vpf);zeros(1,vn.during-vn.vpf)];
    f.f(1,:)=f.f(1,:).*tscale;   %changing row 1 acquisitions to time
    for i=0:tif.xyz(3)-2 %adding appropriate point number in between frames
        f.f(2,(i*vn.vpf)+1:(i+1)*vn.vpf)=linspace(f.fluo(i+1),f.fluo(i+2),vn.vpf);
    end
    
    
    %including timescale to v and isolating only paired recording time
    v.v = [1:length(v.v); v.v];
    v.v(1,:)=v.v(1,:).*tscale;
    v.paired = [f.f(1, :); v.v(2, (vn.pre+vn.vpf+1):(length(v.v)-vn.post))];
    
    
    
    %% outputs
    
    
    %making record struct with recording name
    rec.input_file_names={v.name.name,...
        (strcat('phys samples=',num2str(length(v.v))));...
        tif.name.name,(strcat('image frames=',num2str(tif.xyz(3))));...
        tif.name.name,(strcat('xy pixels=',num2str(tif.xyz(1)),'x,',...
        num2str(tif.xyz(2)),'y'))};
    rec.v=v.v;
    rec.roi_mask=roi.mask;
    rec.background_mask=bk.mask;
    rec.parameters=paras;
    rec.fullframe_reads=tif.reads;
    rec.roi_reads=roi.fluo;
    rec.background_reads=bk.fluo;
    rec.timescale=tscale;
    rec.f=f.f;
    rec.vp=v.paired;
    
    
    %figs
    plot(f.f(1, :), f.f(2, :),'Color',[0.9 0 0])
    %axis([0 1000 23 24.5])
    box off
    %title('Background subtracted fluorescence', 'FontSize',14)
    ax=gca;
    ax.FontSize = 14;
    xlabel(figs.xtitle, 'FontSize',18)
    ylabel('fluorescence', 'FontSize',18)
    savefig(fullfile(path.figs,strcat(path.name,'_fluorescence.fig')))
    
    
    set(figure,'defaultAxesColorOrder',[figs.left_color; figs.right_color],...
        'Color',figs.bk_color);
    yyaxis right
    plot(f.f(1, :), f.f(2, :),'Color',[0.85 0 0])
    %axis([0 100000 -5 15])
    title(path.name, 'FontSize',14)
    ax1=gca;
    ax1.FontSize = 14;
    ylabel('Fluorescence', 'FontSize',12)
    hold on
    yyaxis left
    plot(v.paired(1, :), v.paired(2, :))
    box off
    ax2=gca;
    ax2.FontSize = 14;
    ylabel('mV', 'FontSize',18)
    xlabel(figs.xtitle, 'FontSize',18)
    %axis([0 100000 -5 15])
    hold off
    savefig(fullfile(path.figs,strcat(path.name,'_paired.fig')))
    
    
    close all
    
    
    %saving workspace%
    save(fullfile(path.vars,strcat(path.name,'_workspace.mat')));
    
    
    
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