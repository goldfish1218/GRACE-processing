function example_func_iter50_comparison(input_choice,plot_choice)
%% example_func_iter50_comparison(input_choice,plot_choice)
%% input_choice: at least the first 4 letters of 'constrained' or 'unconstrained'
%% plot_choice: 1,True Model and MTru; 2,M_O_b_s and MPre.
%%
clc
%% prepare input
if strncmp(input_choice,'constrained',4)
    load('constrainedFM.mat','EWHtrue_k_c','EWHpred_k_c')
    
else
    if strncmp(input_choice,'unconstrained',4)
        load('unconstrainedFM.mat','EWHtrue_k_c','EWHpred_k_c')
        else
            error('input_choice cannot be recognized.')
    end
end
pred = cell2mat(EWHpred_k_c');
tru = cell2mat(EWHtrue_k_c');
load('input_var_for_example_Yangtze.mat', 'EWHpred_yzt_SH60_G300', ...
    'model_yzt', ...
    'Lat_Data', ...
    'Lon_Data')
landareas = shaperead('landareas.shp','UseGeoCoords',true);
%%
switch plot_choice
    case 1
        inplot1 = model_yzt;
        inplot2 = reshape(tru(:,50),180,360);
    case 2
        inplot1 = reshape(EWHpred_yzt_SH60_G300,180,360);
        inplot2 = reshape(pred(:,50),180,360);
end

figure,
set(gcf,'units','cent','position',[2,2,24,12]);

subsubplot(1,2,1)
axesm('eqdcylin','frame','on','flinewidth',1,'MapParallels',10,'mlinelocation',10,'grid','on','maplatlim',[15,50],'maplonlim',[80,130])
axis off;
hold on;

geoshow(Lat_Data,Lon_Data,inplot1,'displaytype','texture')
geoshow(landareas,'facecolor','none')
geoshow(Lat_Data(model_yzt==1),Lon_Data(model_yzt==1),'displaytype','point','markeredgecolor','k')

caxis([-1,1])
cmp = colormap('jet');
Ncmp = 100;
cmpnew = [...
    interp1(1:length(cmp(:,1)),cmp(:,1),linspace(1,length(cmp(:,1)),Ncmp))',...
    interp1(1:length(cmp(:,2)),cmp(:,2),linspace(1,length(cmp(:,1)),Ncmp))',...
    interp1(1:length(cmp(:,3)),cmp(:,3),linspace(1,length(cmp(:,1)),Ncmp))'];
colormap(cmpnew)

switch plot_choice
    case 1
        title('(a) true model')
    case 2
        title('(a) M_O_b_s')
end


subsubplot(1,2,2)
axesm('eqdcylin','frame','on','flinewidth',1,'MapParallels',10,'mlinelocation',10,'grid','on','maplatlim',[15,50],'maplonlim',[80,130])
axis off;
hold on;

geoshow(Lat_Data,Lon_Data,inplot2,'displaytype','texture')
geoshow(landareas,'facecolor','none')
geoshow(Lat_Data(model_yzt==1),Lon_Data(model_yzt==1),'displaytype','point','markeredgecolor','k')

caxis([-1,1])
colorbar
colormap(cmpnew)
hcb = colorbar('Position',...
    [0.420576037669742 0.191255787037037 0.195399078070999 0.0470370370370372],'orientation','horiz');
ylabel(hcb,'mEWT')

switch plot_choice
    case 1
        title('(b) M_T_r_u for No.50 iteration')
    case 2
        title('(b) M_P_r_e for No.50 iteration')
end
