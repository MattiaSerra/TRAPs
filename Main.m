%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Author: Mattia Serra                  %
                    % Email:  serram@seas.haravrd.edu       %
                    % Date:   07/07/2020                    %
                    % Web:  https://www.mattiaserra.com/    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The present code computes Transient Attracting Profiles TRAPs [1] [or Attracting Objective Eulerian Coherent Structures(OECSs)] using the theoretical results derived in [2].

% Reference:
%[1] M. Serra, P. Sathe, I. Rypina, A. Kirincich, S. Ross, P. Lermusiaux, A. Allen, T. Peacock and G. Haller, 
% "Search and rescue at sea aided by hidden flow structures",  Nature Communications, 11 2525 (2020).
%[2] M. Serra, and G. Haller, Objective Eulerian Coherent Structures,  
%    Chaos, 26 (2016) 053110.

% Acknowledgements:
% Alireza Hajighasem, Mohammad Farazmand and George Haller.
%%
clear all; close all; clc

% add path with subsuctions and data 
fp = mfilename('fullpath');
    rootdir = fileparts(fp);
    p1 = fullfile(rootdir,'src');
    addpath(rootdir,p1);
    p2 = fullfile(rootdir,'data');
    addpath(rootdir,p2);

% Inputs 
mode = 'serial';                                %'serial' if using a single core, 'parallel' if using multiple cores
options = odeset('RelTol',1e-4,'AbsTol',1e-4);  % ODE solver options
arclength = 1;                                  % Upper bound for TRAP length
NumPointsOnCurve = 40;                          % Number of points along each TRAP
FracAttr = .3;                                  % Requires that the TRAP is attracting at least FracAttr*(atrtraction at the core)


% load velocity field data  
load DataTRAPs.mat 
    % ug:  u velocity
    % vg:  v velocity
    % XG:  x grid 
    % YG:  y grid 

% COMMENT1: This dataset is for illustration. It correspods to a mesoscale
% ocean velocity field produced by SSALTO/DUACS and distributed by AVISO, with support from CNES (http://www.aviso.oceanobs.com/duacs). 

% Compute the eigenvalues and eigenvectors of the rate of stain tensor 
[s,e] = ComputeS_s_e(XG,YG,ug,vg);
    % s: eigenvalues of S.  The first entry denotes the eigenvalue number, the last two entries the grid position 
    % e: eigenvectors of S. The first entry denotes the eigenvector number, the second entry denotes the component (x=1, y=2),
                          % the last two entries the grid position

% COMMENT2: Taking spatial derivatives of the velocity field must be done
% accounting for the way the velocity field is computed. For exampe, if
% using High-Frequency-Radar velocity fields, velocity variations within a
% minimum length scale cannot be resolved. This tipically follows from the algorithm used to reconstract the velocity field from HFR signals (see [2] for details). 

% Smallest eigenvalue 
s1 = squeeze(s(1,:,:));
%% find extema of s1 
[XXtot,YYtot,flag]=findSurfextrema(XG,YG,s1);
s1tot = interp2(XG,YG,s1,XXtot,YYtot,'spline',0);
    % XXtot x coordinate of the critical points 
    % YYtot y coordinate of the critical points 
    % flag = 1; % minima flag = 2; % saddle flag = 3; % maxima
    % s1tot = s1 value at the critical points 

% COMMENT3: With noisy velocity field, using the Hessian method to identify minima of s1 can be inaccurate. 
% In these cases, please use more appropriate algorithms. 

% TRAP Cores 
xTC = XXtot(flag==1);  % x coordinate of the TRAP cores 
yTC = YYtot(flag==1);  % y coordinate of the TRAP cores 
s1TC = s1tot(flag==1); % normal attraction rate at the TRAP cores 

        % Plot s1 contours 
        hcont = figure('units','normalized','outerposition',[0 0 1 1]);
        contour(XG(1,:),YG(:,1),s1,100);shading interp; hold on
        xlabel('$$x$$','Interpreter','latex','FontWeight','bold','FontSize',20);ylabel('$$y$$','Interpreter','latex','FontWeight','bold','FontSize',20);
        axis equal tight; set(gcf,'color','w'); set(gca,'FontSize',16,'fontWeight','normal');set(gca,'YDir','normal')
        title('TRAPs and attraction rate map: s_{1}','Interpreter','tex','FontWeight','bold','FontSize',20)
        colormap(gca,'parula'); hhF=colorbar(gca); hhF.Location='eastoutside'; hhF.FontSize=20;
        set(get(hhF,'ylabel'),'string','$$s_{1}$$','Interpreter','latex','FontWeight','bold');
        % Plot TRAPs core 
        for kkCP=1:length(xTC)
           plot(xTC(kkCP),yTC(kkCP),'Color','r','Marker','o','Markersize',6,'MarkerFaceColor','r','MarkerEdgeColor','r')  
        end

%% Compute TRAPs
% Define the arclength vector for each TRAP
ArcLength = linspace(0,arclength/2,NumPointsOnCurve/2);

% Compute tensorlines of the dominant eigenvector field of the rate of
% strain tensor starting at minima of s1
[pxt_raw,pyt_raw] = TRAPs(xTC,yTC,ArcLength,XG,YG,e,mode,options);

% Make sure that TRAPs are attracting and their normal attraction is greater than a desired fraction of the attraction at the core
[pxt,pyt] = TRAPs_process(pxt_raw,pyt_raw,XG,YG,s1,FracAttr);

        % Plot Traps 
        for kk = 1:size(pxt,2)
        plot(pxt(:,kk),pyt(:,kk),'Color','r','Linewidth',2)
        end 
