%%% visContactMap.m
%%% visulize the simulated contact map vs. Hi-C map 
%%% Last modified: 04/24/2018

%% information about the input parameters:
% celltype            % options refer to README 
% chrId               % chromosome id, options refer to README
% resolution          % resolution to visualize contact map
                      % Hi-C map should have consistent resolution
% startnum            % starting postion of the genomic segment
% endnum              % ending postion of the genomic segment
% nbead               % number of beads included to visualize
                      % the size of the contact matrix
% simmap_path         % path to simulated contact map
% hic_path            % path to raw Hi-C map, 
                      % Rao, Suhas S.P. et al. Cell 159, 1665-1680 (2014).
% hic_norm_path       % path to normalization constant of Hi-C map

%% clear all
clc;clf;clear;close all
addpath('../../src/cmap/MATLAB/')

%% initialize color map (red, green, blue)
load MyColormaps cmap64;
cmap256 = zeros(1024,3);
cmap256(:,1) = interp1(linspace(0,1,64),cmap64(:,1), linspace(0,1,1024));
cmap256(:,2) = interp1(linspace(0,1,64),cmap64(:,2), linspace(0,1,1024));
cmap256(:,3) = interp1(linspace(0,1,64),cmap64(:,3), linspace(0,1,1024));
set(gcf, 'Colormap', cmap256)

%% Initialization
% initialize the param for sim contact map
fprintf('\n**** Processing simulated map ****\n\n')
[simmap,lbl,lbl_tick,...
 celltype,chrId,resolution,nbead,...
 hic_path,hic_norm_path,...
 gpSta,gpEnd] = initcmap(param);

%% Hi-C map
% initialize the param for hic map
hicname = ['./hic/hicMat/' celltype '_chr' num2str(chrId) ...
           '_' num2str(resolution/1E3) 'kb_' num2str(gpSta) ...
           '_' num2str(gpEnd) '.mat'];
if ~exist(hicname,'file')
    fprintf('\n**** Processing Hi-C map ****\n\n')
    % process the hic data
    conv2mat(celltype,chrId,resolution)
    balanceHiC(celltype,chrId,resolution,gpSta,gpEnd)
    normConst(celltype,chrId,resolution,gpSta,gpEnd)
end

%% combine Hi-C with simulated map
cmap = combcmap(nbead,hic_path,hic_norm_path,simmap);

%% Simulated map
% plot the simulated contact map
imagesc(cmap)
colormap(cmap256)
caxis([-7.5,0.1])
colorbar

%% image settings
set(gcf,'Position',[0 0 600 500])
xticks(lbl);
xticklabels(lbl_tick);
yticks(lbl);
yticklabels(lbl_tick);
title([celltype,', chromosome ' num2str(chrId)])
xlabel ('Genomic sequence (Mb)','FontWeight','bold');
ylabel ('Genomic sequence (Mb)','FontWeight','bold');
set(gca,'FontSize',26)
