function example_func_iteration_comparison(input_choice,plot_choice)
%% iteration_comparison(input_choice,plot_choice)
%% input_choice: at least the first 4 letters of 'constrained' or 'unconstrained'
%% plot_choice: 1,M_Pre minus MObs; 2,True Model minus MTru;3,M_Pre;4,MTru.
%%
clc
%% prepare input
if strncmp(input_choice,'constrained',4)
    load('constrainedFM.mat','EWHtrue_k_c','EWHpred_k_c','RMSE_k')
else
    if strncmp(input_choice,'unconstrained',4)
        load('unconstrainedFM.mat','EWHtrue_k_c','EWHpred_k_c','RMSE_k')
        else
            error('input_choice cannot be recognized.')
    end
end
%% plot
landareas = shaperead('landareas.shp','UseGeoCoords',true);
load('input_var_for_example_Yangtze.mat', 'EWHpred_yzt_SH60_G300', ...
    'model_yzt', ...
    'Lat_Data', ...
    'Lon_Data')
figure,
set(gcf,'units','cent','position',[2,2,24,20]);
pred = cell2mat(EWHpred_k_c');
tru = cell2mat(EWHtrue_k_c');
for iii = 1:9
    switch plot_choice
        case 1
            inplot = reshape(EWHpred_yzt_SH60_G300-pred(:,iii),180,360);
            cax_lim = 0.2;
        case 2
            inplot = model_yzt-reshape(tru(:,iii),180,360);
            cax_lim = 0.5;
        case 3
            inplot = reshape(pred(:,iii),180,360);
            cax_lim = 1;
        case 4
            inplot = reshape(tru(:,iii),180,360);
            cax_lim = 1;
    end
    subsubplot(3,3,iii,'vpad',0.01,'hpad',0.01)
    axesm('eqdcylin','frame','on','flinewidth',1,'MapParallels',10,'mlinelocation',10,'grid','on','maplatlim',[15,50],'maplonlim',[80,130])
    axis off;
    hold on;
    
    geoshow(Lat_Data,Lon_Data,inplot,'displaytype','texture')
    geoshow(landareas,'facecolor','none')
    
    caxis([-1,1]*cax_lim)
    cmp = colormap('jet');
    Ncmp = 100;
    cmpnew = [...
        interp1(1:length(cmp(:,1)),cmp(:,1),linspace(1,length(cmp(:,1)),Ncmp))',...
        interp1(1:length(cmp(:,2)),cmp(:,2),linspace(1,length(cmp(:,1)),Ncmp))',...
        interp1(1:length(cmp(:,3)),cmp(:,3),linspace(1,length(cmp(:,1)),Ncmp))'];
    colormap(cmpnew)
    title(['iteration:',num2str(iii),' RMSE:',num2str(RMSE_k(iii))])
end
switch plot_choice
    case 1
        sgtitle('M_P_r_e minus M_O_b_s')
    case 2
        sgtitle('True Model minus M_T_r_u')
    case 3
        sgtitle('M_P_r_e')
    case 4
        sgtitle('M_T_r_u')
end

hcb = colorbar('Position',...
    [0.420576037669742 0.0859236936746601 0.198383541513593 0.0190476189766611],'orientation','horiz');
ylabel(hcb,'mEWT')
