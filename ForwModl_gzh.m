function [EWHtrue_k_c,EWHpred_k_c,RMSE_k] = ForwModl_gzh(EWHobs,EWHtrue_ini,Lat_Data,Lon_Data,land_mask,grdarea,target_mask,RMSE_threshold,iteration_threshold,GauR,maxSHod)
%% [EWHtrue_k_c,EWHpred_k_c,RMSE_k] = ForwModl_gzh(EWHobs,EWHtrue_ini,Lat_Data,Lon_Data,land_mask,grdarea,target_mask,RMSE_threshold,iteration_threshold,GauR,maxSHod)
%% 
%% #input:
%%  EWHobs, observation in unit of m.
%%  EWHtrue_ini, initial model in unit of m.
%%  Lat_Data, latitude in deg.
%%  Lon_Data, longitude in deg.
%%  land_mask, 1 for land, 0 for ocean.
%%  grdarea, grid area in unit of m^2.
%%  target_mask, target region with 1 for target, 0 for outside.
%%  RMSE_threshold, stop criterion for RMSE,in unit of Gt
%%  iteration_threshold, maximum iteration number.
%%  GauR, radius for gaussian filter.
%%  maxSHod, maximum SH degree and order.
%% 
%% #output:
%%  EWHtrue_k_c, modeled true EWH for all iterations,
%%  EWHpred_k_c, modeled prediction for all iterations,
%%  RMSE_k, RMSE for all iterations,
clc
%% #####################prepare##################### 
EWHtrue_k = EWHtrue_ini(:);
land_ind = (land_mask==1);
ocean_ind = (land_mask==0);
target_ind = find(target_mask==1);
land_noObjt_ind = (land_mask==1&target_mask==0);

areaocean = sum(grdarea(ocean_ind));
Mtrue_land = EWHtrue_k(land_ind)'*grdarea(land_ind);
EWHtrue_k(ocean_ind) = -Mtrue_land/areaocean;
RMSE = inf;
Npoint = length(target_ind);
iteration = 0;
EWHtrue_k_c = cell(iteration_threshold+1,1);
EWHtrue_k_c{1} = EWHtrue_k; % ocean=-land
EWHpred_k_c = cell(iteration_threshold,1);
RMSE_k = zeros(iteration_threshold,1);
nSH = (maxSHod+2)*(maxSHod+1)/2;

gain_factor = 1.2;

disp('computing G_SH2EWT')
[G_SH2EWT,cJ_g3d,sJ_g3d,Pnm3d] = inner_make_G(GauR,maxSHod,2,Lon_Data,Lat_Data);
disp('computing G_EWT2SH')
G_EWT2SH = inner_G_grid2SH(Lat_Data,Lon_Data,grdarea*1e-6,maxSHod,Pnm3d,cJ_g3d,sJ_g3d);

reverseStr = '';
%% #####################iteration#####################
while RMSE>RMSE_threshold
    if iteration >iteration_threshold
        fprintf('\n')
        disp(['iteration has reached maximum: ',num2str(iteration_threshold),'. '])
        break;
    end
    iteration = iteration+1;
    % #####################SH transform#####################
    SH_k = G_EWT2SH*EWHtrue_k(:);
    tmp = reshape(SH_k,nSH,2);
    SH_k_C = tmp(:,1);
    SH_k_S = tmp(:,2);
    SH_k_C(1) = 0; % C00=0
    % #####################SH2EWT#####################
    EWHpred_k = G_SH2EWT*[SH_k_C;SH_k_S];
    EWHpred_k_c{iteration} = EWHpred_k; % before correct
    % #####################adjustment#####################
    EWHdiff_k = EWHobs(:)-EWHpred_k;
    EWHtrue_k(target_ind) = EWHtrue_k(target_ind)+gain_factor*EWHdiff_k(target_ind);
    EWHtrue_k(land_noObjt_ind) = 0;
    % #####################update#####################
    EWHtrue_k_c{iteration+1} = EWHtrue_k; % after correct
    % #####################compute RMSE#####################
    RMSE = sqrt(sum((EWHdiff_k(target_ind).*grdarea(target_ind)*1e-9).^2)/Npoint*Npoint);
    RMSE_k(iteration) = RMSE;

    if iteration>1
        msg = sprintf('Processed No.%d/%d RMSE:%.4E', iteration, iteration_threshold,RMSE);
        fprintf([reverseStr,msg]);
        reverseStr = repmat(sprintf('\\b'), 1, length(msg));
    end
    
end
fprintf('\n Computation completed. \n')
end % <main functions> end

%% #####################sub functions#####################
%% List:
%%  G_grid2SH
%%  cal_csGrid
%%  pz_legendre
%%  pz_parpool
%%  make_G
%%  cal_GaussianCoe
%%  cal_EWTxs
%%  love_num

function [G,poolobj,varargout] = inner_G_grid2SH(Lat_Data,Lon_Data,grdarea,maxSHod,varargin)
%% [G,poolobj,varargout] = inner_G_grid2SH(Lat_Data,Lon_Data,grdarea,maxSHod,varargin)
%%    Pnm3d = varargin{1};
%%    cJ_g3d = varargin{2};
%%    sJ_g3d = varargin{3};
%%
disp('using innter G_grid2SH')
rho_ave = 5517;
rho_wat = 1000;
NSH = (1+maxSHod+1)*(maxSHod+1)/2;
ER = 6.3781364600E+06; % variable: mean equator radius; unit: m
if nargout>3 % nargout>3
    extraoutflag = 1;
else
    extraoutflag = 0;
end
if ~isempty(varargin)
    extrainflag = 1;
    Pnm3d = varargin{1};
    cJ_g3d = varargin{2};
    sJ_g3d = varargin{3};
else
    extrainflag = 0;
end
%%
r_ind = (tril(ones(maxSHod+1,maxSHod+1))==1);

n = 0:maxSHod;
[~,n_g1] = meshgrid(n,n);
n_g = n_g1;
n_g = tril(n_g);
n_vec = n_g(r_ind); 

if maxSHod<100
    interp_love_v = love_num(100);
else
    error('exceed the maximum')
end
kn = (interp_love_v(n_vec+1));

r_G = grdarea*1.0E6/ER^2;
l_G = 3*(ones(NSH,1)+kn)./(4*pi*ER*rho_ave/rho_wat*(2*n_vec+1));
if size(l_G,1)>1
    l_G = l_G';
end

[Nrow,Ncol] = size(Lat_Data);
Gup = zeros(Nrow,Ncol,NSH);
Glo = zeros(Nrow,Ncol,NSH);
poolobj = pz_parpool;

tic1 = tic;
if extrainflag
    parfor iii = 1:Nrow
    r_G_in = r_G(iii,:);
    q0 = squeeze(Pnm3d(iii,:,:));
    q1 = squeeze(cJ_g3d(iii,:,:));
    q2 = squeeze(sJ_g3d(iii,:,:));
    r_G_in = repmat(r_G_in,[NSH,1]);
    r_G_in = r_G_in';
    Gup(iii,:,:) = l_G.*q1.*q0.*r_G_in;
    Glo(iii,:,:) = l_G.*q2.*q0.*r_G_in;
%     disp(strcat(num2str(round(iii/Nrow*100,1)),'% has been done'))
%     toc(tic1)
    end % parfor end
else
    Pnm3d = zeros(Nrow,Ncol,NSH);
    cJ_g3d = zeros(Nrow,Ncol,NSH);
    sJ_g3d = zeros(Nrow,Ncol,NSH);
    parfor iii = 1:Nrow
        
        r_G_in = r_G(iii,:);
        inlat = Lat_Data(iii,:);

        q0 = pz_legendre(inlat,maxSHod);
        inlon = Lon_Data(iii,:);
        [cJ_g3d_row,sJ_g3d_row] = cal_csGrid(deg2rad(inlon),maxSHod);
        q1 = squeeze(cJ_g3d_row);
        q2 = squeeze(sJ_g3d_row);

        r_G_in = repmat(r_G_in,[NSH,1]);
        r_G_in = r_G_in';
        Gup(iii,:,:) = l_G.*q1.*q0.*r_G_in;
        Glo(iii,:,:) = l_G.*q2.*q0.*r_G_in;
%         toc(tic1)
        if extraoutflag
            Pnm3d(iii,:,:) = q0;
            cJ_g3d(iii,:,:) = cJ_g3d_row;
            sJ_g3d(iii,:,:) = sJ_g3d_row;
        end % extraoutflag end
%         disp(strcat('SH transform',32,num2str(round(iii/Nrow*100,1)),'% has been done'))
    end % parfor end
end % extrainflag end
% delete(poolobj);
toc(tic1)
disp(strcat('SH transform has been done'))
G = cat(3,Gup,Glo);
G = reshape(G,Nrow*Ncol,2*NSH);
G = G';
if extraoutflag
    varargout{1} = Pnm3d;
    varargout{2} = cJ_g3d;
    varargout{3} = sJ_g3d;
end % extraoutflag end
end % <G_grid2SH> end

function [cJ_g3d,sJ_g3d] = cal_csGrid(J_grid,SHod)
% unit of J : rad
%%
r_ind = find(tril(ones(SHod+1,SHod+1)));
NR = size(J_grid,1);
NC = size(J_grid,2);

if all(all(~diff(J_grid,1,1)))
    copyRows = NR;
    NR=1;
    repeat_flag = 1;
else
    if all(all(~diff(J_grid,1,2)))
        copyCols = NC;
        NC=1;
        repeat_flag = 2;
    else
        repeat_flag = 0;
    end
end
% repeat_flag
cJ_g3d = zeros(NR,NC,length(r_ind));
sJ_g3d = zeros(NR,NC,length(r_ind));

m = 0:SHod;
[m_g1,~] = meshgrid(m,m);
m_g = m_g1;
m_g = tril(m_g);
m_vec = m_g(r_ind);

for iii = 1:NR
    for ii = 1:NC
        J = J_grid(iii,ii);
        cJ_g3d(iii,ii,:) = cos(J*m_vec);
        sJ_g3d(iii,ii,:) = sin(J*m_vec);
    end
end
switch repeat_flag
    case 0
    case 1
        cJ_g3d = repmat(cJ_g3d,copyRows,1,1);
        sJ_g3d = repmat(sJ_g3d,copyRows,1,1);
    case 2
        cJ_g3d = repmat(cJ_g3d,1,copyCols,1);
        sJ_g3d = repmat(sJ_g3d,1,copyCols,1);
end
cJ_g3d = squeeze(cJ_g3d);
sJ_g3d = squeeze(sJ_g3d);
end % <cal_csGrid> end

function [Pmat] = pz_legendre(latitude,maxSHod)
%%
if size(latitude,2)>1
    latitude = latitude';
end
W = latitude;
NNlat = length(W);
t = sind(W);
u = cosd(W);
inner = zeros(maxSHod+1,maxSHod+1,NNlat);
inner(1,1,:) = 1;
inner(2,1,:) = sqrt(3)*t;
inner(2,2,:) = sqrt(3)*u;
for m=0:maxSHod
    for n = m:maxSHod
        if m<2
            if m==n && m >= 2
                Pmm = squeeze(inner(m,m,:));
                inner(m+1,m+1,:) =  sqrt((2*m+1)/(2*m))*u.*Pmm;
            else
                if n>1
                    Pnm1 = squeeze(inner(n,m+1,:));
                    Pn1m1 = squeeze(inner(n-1,m+1,:));
                    inner(n+1,m+1,:) = sqrt((2*n-1)*(2*n+1)/((n-m)*(n+m)))*t.*Pnm1-sqrt((2*n+1)*(n+m-1)*(n-m-1)/((2*n-3)*(n+m)*(n-m))).*Pn1m1;
                end
            end
        else
            alpha1 = sqrt((2*n+1)*(n-m)*(n-m-1)/(2*n-3)/(n+m)/(n+m-1));
            if m==2
                gamma1 = sqrt(2)*sqrt((n-m+1)*(n-m+2)/(n+m)/(n+m-1));
                beta1 = sqrt(2)*sqrt((2*n+1)*(n+m-2)*(n+m-3)/(2*n-3)/(n+m)/(n+m-1));
            else
                gamma1 = sqrt((n-m+1)*(n-m+2)/(n+m)/(n+m-1));
                beta1 = sqrt((2*n+1)*(n+m-2)*(n+m-3)/(2*n-3)/(n+m)/(n+m-1));
            end
            Pnm1m1 = squeeze(inner(n-1,m+1,:));
            Pnm1mm1 = squeeze(inner(n-1,m-1,:));
            Pn1m1 = squeeze(inner(n+1,m-1,:));
                inner(n+1,m+1,:) = alpha1.*Pnm1m1+beta1.*Pnm1mm1-gamma1.*Pn1m1;
        end
    end
end
NSH = (1+maxSHod+1)*(maxSHod+1)/2;
Pmat = zeros(length(W),NSH);
r_ind = (tril(ones(maxSHod+1,maxSHod+1))==1);
inner = reshape(inner,(maxSHod+1)^2,NNlat);
for iii = 1:NNlat
    Pmat(iii,:) = inner(r_ind,iii);
end
end % <pz_legendre> end

function poolobj = pz_parpool()
if isempty(gcp('nocreate'))
    ncores = feature('numCores'); %
    poolobj = parpool('local',ncores);
    disp(strcat('matlab pool start and uses',32,num2str(ncores),32,'cores'));
else
    poolobj = gcp('nocreate');
    disp('matlab pool already started');
end
end % <pz_parpool> end

function [G,cJ_g3d,sJ_g3d,Pnm3d] = inner_make_G(Gau_radii,SHdeg,varFlag,varargin)
%% [G,cJ_g3d,sJ_g3d,Pnm3d] = inner_make_G(Gau_radii,SHdeg,varFlag,varargin)
%% case 1
%%        Pnm3d = varargin{1};
%%        cJ_g3d = varargin{2};
%%        sJ_g3d = varargin{3};
%% case 2
%%        Lon_Data = varargin{1};
%%        Lat_Data = varargin{2};
%% ===============================
disp('using innter make_G')
switch varFlag
    case 1
        Pnm3d = varargin{1};
        cJ_g3d = varargin{2};
        sJ_g3d = varargin{3};
        [nrow,ncol,nSH] = size(Pnm3d);
    case 2
        nSH = (1+SHdeg+1)*(SHdeg+1)/2;
        Lat_Data = varargin{2};
        Lon_Data = varargin{1};
        [nrow,ncol] = size(Lat_Data);
        Pnm3d = zeros(nrow,ncol,nSH);
        cJ_g3d = zeros(nrow,ncol,nSH);
        sJ_g3d = zeros(nrow,ncol,nSH);
        for iii = 1:nrow
            inlat = Lat_Data(iii,:);
            q0 = pz_legendre(inlat,SHdeg);
            inlon = Lon_Data(iii,:);
            [cJ_g3d_row,sJ_g3d_row] = cal_csGrid(deg2rad(inlon),SHdeg);
            Pnm3d(iii,:,:) = q0;
            cJ_g3d(iii,:,:) = cJ_g3d_row;
            sJ_g3d(iii,:,:) = sJ_g3d_row;
        end
end
%%
W_filter = cal_GaussianCoe(Gau_radii,100);
Pnm2d = reshape(Pnm3d,nrow*ncol,nSH);
cJ_g2d = reshape(cJ_g3d,nrow*ncol,nSH);
sJ_g2d = reshape(sJ_g3d,nrow*ncol,nSH);
EWT_xsVec = cal_EWTxs(W_filter,SHdeg);
G = [Pnm2d,Pnm2d].*[cJ_g2d,sJ_g2d].*repmat((EWT_xsVec)',[nrow*ncol,2]);
end % <make_G> end

function W_filter = cal_GaussianCoe(R,deg)
%%
ER = 6.3781364600E+06;
Gau_r = R*1.0E3;
Gau_b = log(2)/(1-cos(Gau_r/ER));
W_filter(1) = 1;
W_filter(2) = ((1+exp(-2*Gau_b))/(1-exp(-2*Gau_b))-1/Gau_b);
for i = 3:deg
    W_filter(i) = (-(2*i-3)/Gau_b*W_filter(i-1)+W_filter(i-2));
end
% disp('生成高斯滤波因子')
end % <cal_GaussianCoe> end

function EWT_xsVec = cal_EWTxs(W_filter,SHod)
%%
rho_wat = 1000;
rho_ave = 5517;
ER = 6.3781363000E+06;
r_ind = (tril(ones(SHod+1,SHod+1))==1);
n = 0:SHod;
[~,n_g1] = meshgrid(n,n);
n_g = n_g1;
n_g = tril(n_g);
n_vec = n_g(r_ind);
interp_love_v = love_num(100);
kn = (interp_love_v(n_vec+1));
Wn = (W_filter(n_vec+1));
if size(kn,2)>1
    kn = kn';
end
if size(Wn,2)>1
    Wn = Wn';
end
if size(n_vec,2)>1
    n_vec = n_vec';
end
EWT_xsVec = ((ER*rho_ave/3/rho_wat)*Wn.*((2*n_vec+1)./(1+kn)));
end % <cal_EWTxs> end

function output = love_num(Trun_deg)
love_all = [0;0.0270000000000000;-0.303000000000000;-0.194000000000000;-0.132000000000000;-0.104000000000000;-0.0890000000000000;-0.0810000000000000;-0.0760000000000000;-0.0720000000000000;-0.0690000000000000;-0.0665000000000000;-0.0640000000000000;-0.0620000000000000;-0.0600000000000000;-0.0580000000000000;-0.0566000000000000;-0.0552000000000000;-0.0538000000000000;-0.0524000000000000;-0.0510000000000000;-0.0499000000000000;-0.0488000000000000;-0.0477000000000000;-0.0466000000000000;-0.0455000000000000;-0.0444000000000000;-0.0433000000000000;-0.0422000000000000;-0.0411000000000000;-0.0400000000000000;-0.0393000000000000;-0.0386000000000000;-0.0379000000000000;-0.0372000000000000;-0.0365000000000000;-0.0358000000000000;-0.0351000000000000;-0.0344000000000000;-0.0337000000000000;-0.0330000000000000;-0.0324000000000000;-0.0318000000000000;-0.0312000000000000;-0.0306000000000000;-0.0300000000000000;-0.0294000000000000;-0.0288000000000000;-0.0282000000000000;-0.0276000000000000;-0.0270000000000000;-0.0266500000000000;-0.0263000000000000;-0.0259500000000000;-0.0256000000000000;-0.0252500000000000;-0.0249000000000000;-0.0245500000000000;-0.0242000000000000;-0.0238500000000000;-0.0235000000000000;-0.0231500000000000;-0.0228000000000000;-0.0224500000000000;-0.0221000000000000;-0.0217500000000000;-0.0214000000000000;-0.0210500000000000;-0.0207000000000000;-0.0203500000000000;-0.0200000000000000;-0.0198000000000000;-0.0196000000000000;-0.0194000000000000;-0.0192000000000000;-0.0190000000000000;-0.0188000000000000;-0.0186000000000000;-0.0184000000000000;-0.0182000000000000;-0.0180000000000000;-0.0178000000000000;-0.0176000000000000;-0.0174000000000000;-0.0172000000000000;-0.0170000000000000;-0.0168000000000000;-0.0166000000000000;-0.0164000000000000;-0.0162000000000000;-0.0160000000000000;-0.0158000000000000;-0.0156000000000000;-0.0154000000000000;-0.0152000000000000;-0.0150000000000000;-0.0148000000000000;-0.0146000000000000;-0.0144000000000000;-0.0142000000000000;-0.0140000000000000;-0.0139200000000000;-0.0138400000000000;-0.0137600000000000;-0.0136800000000000;-0.0136000000000000;-0.0135200000000000;-0.0134400000000000;-0.0133600000000000;-0.0132800000000000;-0.0132000000000000;-0.0131200000000000;-0.0130400000000000;-0.0129600000000000;-0.0128800000000000;-0.0128000000000000;-0.0127200000000000;-0.0126400000000000;-0.0125600000000000;-0.0124800000000000;-0.0124000000000000;-0.0123200000000000;-0.0122400000000000;-0.0121600000000000;-0.0120800000000000;-0.0120000000000000;-0.0119200000000000;-0.0118400000000000;-0.0117600000000000;-0.0116800000000000;-0.0116000000000000;-0.0115200000000000;-0.0114400000000000;-0.0113600000000000;-0.0112800000000000;-0.0112000000000000;-0.0111200000000000;-0.0110400000000000;-0.0109600000000000;-0.0108800000000000;-0.0108000000000000;-0.0107200000000000;-0.0106400000000000;-0.0105600000000000;-0.0104800000000000;-0.0104000000000000;-0.0103200000000000;-0.0102400000000000;-0.0101600000000000;-0.0100800000000000;-0.0100000000000000;-0.00994000000000000;-0.00988000000000000;-0.00982000000000000;-0.00976000000000000;-0.00970000000000000;-0.00964000000000000;-0.00958000000000000;-0.00952000000000000;-0.00946000000000000;-0.00940000000000000;-0.00934000000000000;-0.00928000000000000;-0.00922000000000000;-0.00916000000000000;-0.00910000000000000;-0.00904000000000000;-0.00898000000000000;-0.00892000000000000;-0.00886000000000000;-0.00880000000000000;-0.00874000000000000;-0.00868000000000000;-0.00862000000000000;-0.00856000000000000;-0.00850000000000000;-0.00844000000000000;-0.00838000000000000;-0.00832000000000000;-0.00826000000000000;-0.00820000000000000;-0.00814000000000000;-0.00808000000000000;-0.00802000000000000;-0.00796000000000000;-0.00790000000000000;-0.00784000000000000;-0.00778000000000000;-0.00772000000000000;-0.00766000000000000;-0.00760000000000000;-0.00754000000000000;-0.00748000000000000;-0.00742000000000000;-0.00736000000000000;-0.00730000000000000;-0.00724000000000000;-0.00718000000000000;-0.00712000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000];
output = love_all(1:Trun_deg+1);
end % <love_num> end
