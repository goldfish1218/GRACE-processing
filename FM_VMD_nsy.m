function [M_ini_k_c,M_fm_k_c,RMSE_v,m_land_c,RMSE_n] = FM_VMD_nsy(S_obs,N_obs,lat,lon,land_mask,grdarea,object_mask,RMSE_threshold,iteration_threshold,varargin)
%%  [M_ini_k_c,M_fm_k_c,RMSE_v,m_land_c] = FM_VMD_nsy(S_obs,N_obs,lat,lon,land_mask,grdarea,object_mask,RMSE_threshold,iteration_threshold,varargin)
%%  
%%  iteration_threshold = 100;
%%  1,
%%
if ~isempty(varargin)
%         Pnm3d = varargin{1};
%         cJ_g3d = varargin{2};
%         sJ_g3d = varargin{3};
%         G = varargin{4};
    G = varargin{1};
    invG = varargin{2};
%         invG = pinv(G);
    iniflag = 1;
else
    iniflag = 0;
end
M_obs = S_obs(:);
S_ini_k = S_obs(:);

% when the ocean mass in initial model does not balance the land mass
% to conserve the total mass.
% if do.
% then
land_ind = find(land_mask);
ocean_ind = find(land_mask==0);
arealand = sum(grdarea(land_ind));
areaocean = sum(grdarea(ocean_ind));
Sland = S_ini_k(land_ind)'*grdarea(land_ind);
S_ini_k(ocean_ind) = -Sland/areaocean;

% add noise
M_ini_k = S_ini_k;

%     dRMSE = 9999;
RMSE = 9999;
Npoint = length(M_ini_k);
iteration = 0;
M_ini_k_c = cell(iteration_threshold,1);
M_fm_k_c = cell(iteration_threshold,1);
m_land_c = zeros(iteration_threshold,1);
RMSE_v = zeros(iteration_threshold,1);
RMSE_n = zeros(iteration_threshold,1);
[Nrow,Ncol] = size(land_mask);
%%
%%  2,
%%
    gain_factor = 1.5*ones(iteration_threshold+1,1);
% gain_factor = 1.2.^linspace(0,0,101);
object_ind = find(object_mask==1);
sf_noise = 1;
reverseStr = '';
while RMSE>RMSE_threshold
    if iteration >iteration_threshold
        fprintf('\n')
        disp('iteration has reached maximum: 100. might not reach convergence')
        break
    end
    iteration = iteration+1;
    M_ini_k_c{iteration} = M_ini_k;

%         tran_M_ini_k = reshape(M_ini_k,Nrow,Ncol);


    if iniflag==0
        if iteration==1
%                 [M_fmkc,M_fmks,~,Pnm3d,cJ_g3d,sJ_g3d] = cal_grid2SHeasy(lat,lon,tran_M_ini_k,grdarea,60);
            Pnm3d = zeros(181,361,1891);
            pz_parpool();
            tic1 =tic;
            parfor iii = 1:181
                inlat = lat(iii,:);
                q0 = pz_legendre(inlat,60,3,1);
                Pnm3d(iii,:,:) = q0;
                toc(tic1)
            end
            [cJ_g3d,sJ_g3d] = cal_csGrid(deg2rad(lon),60);
            G = make_G(300,60,1,Pnm3d,cJ_g3d,sJ_g3d);
            invG = pinv(G);
            assignin('base','G',G)
            assignin('base','invG',invG)
        end
    end
%             [M_fmkc,M_fmks] = cal_grid2SHeasy(lat,lon,tran_M_ini_k,grdarea,60,Pnm3d,cJ_g3d,sJ_g3d);
    SH = invG*M_ini_k(:);
    tmp = reshape(SH,1891,2);
    M_fmkc = tmp(:,1);
    M_fmks = tmp(:,2);


    M_fmkc(1) = 0; % C00
%         M_fmkc(2) = 0; % C10
%         M_fmkc(62) = 0; %C11
%         M_fmks(62) = 0; %S11

    M_fm_k = G*[M_fmkc;M_fmks];
    in_vmd = reshape(M_fm_k,Nrow,Ncol)+sf_noise*N_obs;
    % filtering with VMD
%     in_vmd = add_noise+reshape(M_fm_k,Nrow,Ncol);
    [M_signal,M_noise] = BVMD_pz(in_vmd,5,300,0.1,1,0,5*10^-5);
    N_d_k = N_obs(:)-M_noise(:);
    
    S_d_k = S_obs(:)-M_signal(:);
    S_ini_k(object_ind) = S_ini_k(object_ind)+gain_factor(iteration)*S_d_k(object_ind);
    M_fm_k_c{iteration} = M_fm_k;
%         M_fm_k = G*SH;
%     M_d_k = M_obs-M_fm_k;
%     M_ini_k(object_ind) = M_ini_k(object_ind)+gain_factor(iteration)*M_d_k(object_ind);
    M_d_k = S_d_k;
    M_ini_k = S_ini_k;
    
    Mland = M_ini_k(land_ind)'*grdarea(land_ind)/arealand;
    m_land_c(iteration) = Mland;
    sf_increase = Mland/Sland*arealand;
    fprintf('Magnitude of Mass has increased by %.6f\n',sf_increase)
    

    RMSE = sqrt(sum(M_d_k.^2)/Npoint);
    RMSE_v(iteration) = RMSE;
    RMSE_n(iteration) = sqrt(sum(N_d_k.^2)/Npoint);

    if iteration>1
        msg = sprintf('Processed No.%d/%d RMSE:%.6f\n', iteration, iteration_threshold,round(RMSE,6));
        fprintf([reverseStr,msg]);
%             reverseStr = repmat(sprintf('\\b'), 1, length(msg));
    end
    if iteration >9
        if RMSE>RMSE_v(iteration-1)
            fprintf('\n')
            disp('The RMSE has increased in this iteration!!!')
%                 break
        else
            if RMSE==RMSE_v(iteration-1)
                fprintf('\n')
                disp('The RMSE has been unchanged!!!')
%                     break
            else
%                 disp('The RMSE decreases...')
            end
        end
    end


end
fprintf('Computation completed. \n')
end

function [C_vec,S_vec,poolobj,varargout] = cal_grid2SHeasy(Lat_Data,Lon_Data,input_EWT,grdarea,maxSHod,varargin)
%% [C_vec,S_vec,poolobj,varargout] = cal_grid2SHeasy(Lat_Data,Lon_Data,input_EWT,grdarea,maxSHod,varargin)
%% if give varargin
%%    Pnm3d = varargin{1};
%%    cJ_g3d = varargin{2};
%%    sJ_g3d = varargin{3};
%% function include:
%%      pz_legendre
%%      cal_csGrid
%% data include:
%%      'Love_num.mat'
%% units:
%%     input:
%%       Lat Lon: deg
%%       input_EWT: mm
%%       grdarea: km^2
%%     output:
%%       the unit of the expansion of spherical harmonics is mm.
%% 1Gt=10^9*10^3kg
%% =10^12kg
%% mass=EWT*area*rho
%% let EWT=1mm, area=1km^2, rho=1000kg/m^2
%% mass =10^-3m*1000kg/m^3*10^6m^2
%% =10^6m*1kg/m^3*m^2
%% =10^6kg
%% =10^-6Gt
%% i.e.
%% 1mmEWT*1km^2area = 10^-6Gt
%%
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
r_ind = find(tril(ones(maxSHod+1,maxSHod+1)));
C_vec = zeros(1,NSH);
S_vec = zeros(1,NSH);
n = 0:maxSHod;
[~,n_g1] = meshgrid(n,n);
n_g = n_g1;
n_g = tril(n_g);
n_vec = n_g(r_ind); %#ok<FNDSB>
% load('Love_num.mat','interp_love_v')
interp_love_v = love_num(100);
kn = (interp_love_v(n_vec+1));
r_G = input_EWT.*grdarea*1.0E6/ER^2;
l_G = 3*(ones(NSH,1)+kn)./(4*pi*ER*rho_ave/rho_wat*(2*n_vec+1));
if size(l_G,1)>1
    l_G = l_G';
end

[Nrow,Ncol] = size(Lat_Data);

% dispflag = 1;
poolobj = pz_parpool;


tic1 = tic;
if extrainflag
    parfor iii = 1:Nrow
%     if iii==1
%         dispflag =1;
%     end
    r_G_in = r_G(iii,:);
    q0 = squeeze(Pnm3d(iii,:,:));
    q1 = squeeze(cJ_g3d(iii,:,:));
    q2 = squeeze(sJ_g3d(iii,:,:));
    r_G_in = repmat(r_G_in,[NSH,1]);
    r_G_in = r_G_in';
    C_vec = C_vec+sum(l_G.*q1.*q0.*r_G_in,1);
    S_vec = S_vec+sum(l_G.*q2.*q0.*r_G_in,1);
%     if rem(round(iii/Nrow),5)==0
%         disp(strcat('SH transform ',32,num2str(round(iii/Nrow*100,1)),'%...'))
%     end
%     toc(tic1)
    end
else
    parfor iii = 1:Nrow
    %     if iii==1
    %         dispflag =1;
    %     end
        r_G_in = r_G(iii,:);
        inlat = Lat_Data(iii,:);

        q0 = pz_legendre(inlat,maxSHod,3,1);
        inlon = Lon_Data(iii,:);
        [cJ_g3d_row,sJ_g3d_row] = cal_csGrid(deg2rad(inlon),maxSHod);
        q1 = squeeze(cJ_g3d_row);
        q2 = squeeze(sJ_g3d_row);

        r_G_in = repmat(r_G_in,[NSH,1]);
        r_G_in = r_G_in';
        C_vec = C_vec+sum(l_G.*q1.*q0.*r_G_in,1);
        S_vec = S_vec+sum(l_G.*q2.*q0.*r_G_in,1);
        toc(tic1)
        if extraoutflag
            Pnm3d(iii,:,:) = q0;
            cJ_g3d(iii,:,:) = cJ_g3d_row;
            sJ_g3d(iii,:,:) = sJ_g3d_row;
        end
%         if rem(round(iii/Nrow),5)==0
%             disp(strcat('SH transform ',32,num2str(round(iii/Nrow*100,1)),'%...'))
%         end
    end
end
% delete(poolobj);
C_vec = C_vec';
S_vec = S_vec';
if extraoutflag
    varargout{1} = Pnm3d;
    varargout{2} = cJ_g3d;
    varargout{3} = sJ_g3d;
end
end

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
end

function [Pmat] = pz_legendre(latitude,maxSHod,method,input_latvec,varargin)
%% the method3 is the best.
%% the input latitude should be degree.
%% =====================================================================
%% the fully normalized associated legendre function computed by 
%% Matlab function legendre, is not consistent with the real value.
%% for column 1st, i.e. when m=0, the value should be divided by sqrt(1/2)
%% for the rest of the columns, the value should multiply a factor of 2.
%%

% [Pmat] = pz_legendre(latitude,maxSHod,method,input_latgrid,varargin)
if input_latvec
    switch method
        case 1
            if size(latitude,2)>1
                latitude = latitude';
            end
            W = latitude;
            t = sind(W);
            u = cosd(W);
            NNlat = length(W);
            inner = zeros(maxSHod+1,maxSHod+1,NNlat);
            inner(1,1,:) = 1;
            inner(2,1,:) = sqrt(3)*t;
            inner(2,2,:) = sqrt(3)*u;
            for m = 0:maxSHod
                for n = 0:maxSHod
                    if m==n && m >= 2
                        Pmm  = squeeze(inner(m,m,:));
                        inner(m+1,m+1,:) =  sqrt((2*m+1)/(2*m))*u.*Pmm;
                    else
                        if n>1
                            Pnm1 = squeeze(inner(n,m+1,:));
                            Pnm1m1 = squeeze(inner(n-1,m+1,:));
                            inner(n+1,m+1,:) = sqrt((2*n-1)*(2*n+1)/((n-m)*(n+m)))*t.*Pnm1-sqrt((2*n+1)*(n+m-1)*(n-m-1)/((2*n-3)*(n+m)*(n-m))).*Pnm1m1;
                        end
                    end
                end
            end
            NSH = (1+maxSHod+1)*(maxSHod+1)/2;
            Pmat = zeros(length(W),NSH);
            r_ind = find(tril(ones(maxSHod+1,maxSHod+1)));
            inner = reshape(inner,(maxSHod+1)^2,NNlat);
            for iii = 1:NNlat
                Pmat(iii,:) = inner(r_ind,iii);
            end
        case 2
        case 3
%             SHo = varargin{1};
%             SHd = varargin{2};
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
            r_ind = find(tril(ones(maxSHod+1,maxSHod+1)));
            inner = reshape(inner,(maxSHod+1)^2,NNlat);
            for iii = 1:NNlat
                Pmat(iii,:) = inner(r_ind,iii);
            end
    end
else
%% method switch start
switch method
    case 1
%% method 1 . Standard forward column recursion method (Colombo 1981)
%% note that: the zenith angle range from 0 to 180 deg. while the latitude range from -90 to 90
%% the relation is : 
%%                  lat: -90 ~ 90   
%%                  zen: 180 ~ 0
%%               zen = 90 - lat(theta)
%% thus: t = cosd(90-theta) = sind(theta)
%%       u = sind(90-theta) = cosd(theta)
W = latitude;
t = sind(W);
u = cosd(W);
Pmat = zeros(maxSHod+1,maxSHod+1);
Pmat(1,1) = 1;
Pmat(2,1) = sqrt(3)*t;
Pmat(2,2) = sqrt(3)*u;
for m = 0:maxSHod
    for n = 0:maxSHod
        if m==n && m >= 2
            Pmat(m+1,m+1) =  sqrt((2*m+1)/(2*m))*u*Pmat(m,m);
        else
            if n>1
            Pmat(n+1,m+1) = sqrt((2*n-1)*(2*n+1)/((n-m)*(n+m)))*t*Pmat(n,m+1)-sqrt((2*n+1)*(n+m-1)*(n-m-1)/((2*n-3)*(n+m)*(n-m)))*Pmat(n-1,m+1);
            end
        end
    end
end
%% ==============
    case 2
%% method 2. Belikov
W = latitude;
t = sind(W);
u = cosd(W);
% Pmat = zeros(SHod+1,SHod+1);
Pnew = zeros(maxSHod+1,maxSHod+2);
Pnew(1,1) = 1;
for m=0:maxSHod
    for n = m:maxSHod
        if m ==0
            if n>=1
                Pnew(n+1,1) = t*Pnew(n,1)-u/2*Pnew(n,2);
            end
        else
%             if n>=1
%             if n<=(SHod+1)
                Pnew(n+1,m+1) = t*Pnew(n,m+1)-u*(Pnew(n,m+2)/4-Pnew(n,m));
%             else
%                 Pnew(n+1,m+1) = t*Pnew(n,m+1)-u*(Pnew(n,m+2)/4-Pnew(n,m));
%             end
%             end
        end
    end
end
factor1 = zeros(maxSHod+1,maxSHod+1);
factor1(1,1) = 1;
factor1(2,1) = 1;
factor1(2,2) = 1;
for m = 0 :maxSHod
    for n = m:maxSHod
        if n>=2 
            if n-1>=m
                factor1(n+1,m+1) = sqrt(1-m^2/n^2)*factor1(n,m+1);
            else
                factor1(n+1,n+1) = sqrt(1-1/2/n)*factor1(n,n);
            end
        end
    end
end
[~,n_g] = meshgrid(0:maxSHod,0:maxSHod);
factor2 = tril(sqrt(2*n_g+1));
% factor1
% factor2
Pmat = Pnew(1:end,1:end-1).*factor1.*factor2; % .*factor .*factor1.*factor2
%% method 3
%% ===============
    case 3
W = latitude;
t = sind(W);
u = cosd(W);
Pmat = zeros(maxSHod+1,maxSHod+1);
Pmat(1,1) = 1;
Pmat(2,1) = sqrt(3)*t;
Pmat(2,2) = sqrt(3)*u;
for m=0:maxSHod
    for n = m:maxSHod
        if m<2
            if m==n && m >= 2
                Pmat(m+1,m+1) =  sqrt((2*m+1)/(2*m))*u*Pmat(m,m);
            else
                if n>1
                    Pmat(n+1,m+1) = sqrt((2*n-1)*(2*n+1)/((n-m)*(n+m)))*t*Pmat(n,m+1)-sqrt((2*n+1)*(n+m-1)*(n-m-1)/((2*n-3)*(n+m)*(n-m)))*Pmat(n-1,m+1);
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
                Pmat(n+1,m+1) = alpha1*Pmat(n-1,m+1)+beta1*Pmat(n-1,m-1)-gamma1*Pmat(n+1,m-1);
        end
    end
end
                
        
%% method switch end
end
%% end of input_latvec
end
end

function poolobj = pz_parpool()
if isempty(gcp('nocreate'))
    ncores = feature('numCores'); %
    if ncores>6
        ncores=24;
    end
    poolobj = parpool('local',ncores);
    disp(strcat('matlab pool start and uses',32,num2str(ncores),32,'cores'));
else
    poolobj = gcp('nocreate');
    disp('matlab pool already started');
end
end

function G = make_G(Gau_radii,SHdeg,varFlag,varargin)
%% G = make_G(Gau_radii,SHdeg,varFlag,varargin)
%% case 1
%%        Pnm3d = varargin{1};
%%        cJ_g3d = varargin{2};
%%        sJ_g3d = varargin{3};
%% case 2
%%        Lon_Data = varargin{1};
%%        Lat_Data = varargin{2};
%% ===============================
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
            q0 = pz_legendre(inlat,SHdeg,3,1);
            inlon = Lon_Data(iii,:);
            [cJ_g3d_row,sJ_g3d_row] = cal_csGrid(deg2rad(inlon),SHdeg);
            Pnm3d(iii,:,:) = q0;
            cJ_g3d(iii,:,:) = cJ_g3d_row;
            sJ_g3d(iii,:,:) = sJ_g3d_row;
        end
end
%%
W_filter = cal_GaussianCoe(Gau_radii,100);
Pnm3d = reshape(Pnm3d,nrow*ncol,nSH);
cJ_g3d = reshape(cJ_g3d,nrow*ncol,nSH);
sJ_g3d = reshape(sJ_g3d,nrow*ncol,nSH);

EWT_xsVec = cal_EWTxs(W_filter,SHdeg);

G = [Pnm3d,Pnm3d].*[cJ_g3d,sJ_g3d].*repmat((EWT_xsVec)',[nrow*ncol,2]);
end

function W_filter = cal_GaussianCoe(R,deg)
%%
% gaussian smoothing factors in SH field. up to 200 degree
% 出了一个大错，公式使用出错，递推公式中是（2l+1），2，3...阶，对应的是，3，5...，
% 而我写成了（2i-1），2，3...阶，对应的i是，3，4...，对应的是，5，7
% 导致了极不稳定
ER = 6.3781364600E+06;
Gau_r = R*1.0E3;
Gau_b = log(2)/(1-cos(Gau_r/ER));
W_filter(1) = 1;
W_filter(2) = ((1+exp(-2*Gau_b))/(1-exp(-2*Gau_b))-1/Gau_b);
for i = 3:deg
    W_filter(i) = (-(2*i-3)/Gau_b*W_filter(i-1)+W_filter(i-2));
end
disp('生成高斯滤波因子')
end

function EWT_xsVec = cal_EWTxs(W_filter,SHod)
%%
rho_wat = 1000;
rho_ave = 5517;
ER = 6.3781363000E+06;
r_ind = find(tril(ones(SHod+1,SHod+1)));

n = 0:SHod;
[~,n_g1] = meshgrid(n,n);
n_g = n_g1;
n_g = tril(n_g);
n_vec = n_g(r_ind); %#ok<FNDSB>
% load('Love_num.mat','interp_love_v')
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
% [kn,Wn,n_vec] = prepareSurfaceData(kn,Wn,n_vec);
EWT_xsVec = ((ER*rho_ave/3/rho_wat)*Wn.*((2*n_vec+1)./(1+kn)));
end

function output = love_num(Trun_deg)
love_all = [0;0.0270000000000000;-0.303000000000000;-0.194000000000000;-0.132000000000000;-0.104000000000000;-0.0890000000000000;-0.0810000000000000;-0.0760000000000000;-0.0720000000000000;-0.0690000000000000;-0.0665000000000000;-0.0640000000000000;-0.0620000000000000;-0.0600000000000000;-0.0580000000000000;-0.0566000000000000;-0.0552000000000000;-0.0538000000000000;-0.0524000000000000;-0.0510000000000000;-0.0499000000000000;-0.0488000000000000;-0.0477000000000000;-0.0466000000000000;-0.0455000000000000;-0.0444000000000000;-0.0433000000000000;-0.0422000000000000;-0.0411000000000000;-0.0400000000000000;-0.0393000000000000;-0.0386000000000000;-0.0379000000000000;-0.0372000000000000;-0.0365000000000000;-0.0358000000000000;-0.0351000000000000;-0.0344000000000000;-0.0337000000000000;-0.0330000000000000;-0.0324000000000000;-0.0318000000000000;-0.0312000000000000;-0.0306000000000000;-0.0300000000000000;-0.0294000000000000;-0.0288000000000000;-0.0282000000000000;-0.0276000000000000;-0.0270000000000000;-0.0266500000000000;-0.0263000000000000;-0.0259500000000000;-0.0256000000000000;-0.0252500000000000;-0.0249000000000000;-0.0245500000000000;-0.0242000000000000;-0.0238500000000000;-0.0235000000000000;-0.0231500000000000;-0.0228000000000000;-0.0224500000000000;-0.0221000000000000;-0.0217500000000000;-0.0214000000000000;-0.0210500000000000;-0.0207000000000000;-0.0203500000000000;-0.0200000000000000;-0.0198000000000000;-0.0196000000000000;-0.0194000000000000;-0.0192000000000000;-0.0190000000000000;-0.0188000000000000;-0.0186000000000000;-0.0184000000000000;-0.0182000000000000;-0.0180000000000000;-0.0178000000000000;-0.0176000000000000;-0.0174000000000000;-0.0172000000000000;-0.0170000000000000;-0.0168000000000000;-0.0166000000000000;-0.0164000000000000;-0.0162000000000000;-0.0160000000000000;-0.0158000000000000;-0.0156000000000000;-0.0154000000000000;-0.0152000000000000;-0.0150000000000000;-0.0148000000000000;-0.0146000000000000;-0.0144000000000000;-0.0142000000000000;-0.0140000000000000;-0.0139200000000000;-0.0138400000000000;-0.0137600000000000;-0.0136800000000000;-0.0136000000000000;-0.0135200000000000;-0.0134400000000000;-0.0133600000000000;-0.0132800000000000;-0.0132000000000000;-0.0131200000000000;-0.0130400000000000;-0.0129600000000000;-0.0128800000000000;-0.0128000000000000;-0.0127200000000000;-0.0126400000000000;-0.0125600000000000;-0.0124800000000000;-0.0124000000000000;-0.0123200000000000;-0.0122400000000000;-0.0121600000000000;-0.0120800000000000;-0.0120000000000000;-0.0119200000000000;-0.0118400000000000;-0.0117600000000000;-0.0116800000000000;-0.0116000000000000;-0.0115200000000000;-0.0114400000000000;-0.0113600000000000;-0.0112800000000000;-0.0112000000000000;-0.0111200000000000;-0.0110400000000000;-0.0109600000000000;-0.0108800000000000;-0.0108000000000000;-0.0107200000000000;-0.0106400000000000;-0.0105600000000000;-0.0104800000000000;-0.0104000000000000;-0.0103200000000000;-0.0102400000000000;-0.0101600000000000;-0.0100800000000000;-0.0100000000000000;-0.00994000000000000;-0.00988000000000000;-0.00982000000000000;-0.00976000000000000;-0.00970000000000000;-0.00964000000000000;-0.00958000000000000;-0.00952000000000000;-0.00946000000000000;-0.00940000000000000;-0.00934000000000000;-0.00928000000000000;-0.00922000000000000;-0.00916000000000000;-0.00910000000000000;-0.00904000000000000;-0.00898000000000000;-0.00892000000000000;-0.00886000000000000;-0.00880000000000000;-0.00874000000000000;-0.00868000000000000;-0.00862000000000000;-0.00856000000000000;-0.00850000000000000;-0.00844000000000000;-0.00838000000000000;-0.00832000000000000;-0.00826000000000000;-0.00820000000000000;-0.00814000000000000;-0.00808000000000000;-0.00802000000000000;-0.00796000000000000;-0.00790000000000000;-0.00784000000000000;-0.00778000000000000;-0.00772000000000000;-0.00766000000000000;-0.00760000000000000;-0.00754000000000000;-0.00748000000000000;-0.00742000000000000;-0.00736000000000000;-0.00730000000000000;-0.00724000000000000;-0.00718000000000000;-0.00712000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000;-0.00706000000000000];
output = love_all(1:Trun_deg+1);

end

function [M_signal,M_noise] = BVMD_pz(input,K,alpha,tau,DC,init,tol)
%%
% alpha = 300;       % bandwidth constraint
% tau = 0.1;         % Lagrangian multipliers dual ascent time step
% K = 5;              % number of modes
% DC = 1;             % includes DC part (first mode at DC)
% init = 0;           % initialize omegas randomly, may need multiple runs!
% tol = K*10^-5;      % tolerance (for convergence)
%%
Nraw = size(input,1);
Ncol = size(input,2);
% if Ncol>361
%     alpha = Ncol;
% end
fprintf('alpha=%d,tau=%.2f\n',alpha,tau);
[u] = VMD_2D(input, alpha, tau, K, DC, init, tol);
[~,ifSignal] = FFTA4VMD(u);
if all(ifSignal)
    stopf = all(ifSignal);
    iter = 0;
    while stopf==1 && iter<5
        disp('Bad Decomposition, noise cannot be islolated!')
        iter = iter+1;
        fprintf('parameters changed %d/5: alpha=%d,tau=%.2f\n',iter,alpha-iter*50,tau);
        [u] = VMD_2D(input, alpha-iter*50, tau, K, DC, init, tol);
        [~,ifSignal] = FFTA4VMD(u);
        disp(strcat('repeat',32,num2str(iter),32,'times'))
        stopf = all(ifSignal);
    end
    stopf = all(ifSignal);
    if stopf
        disp('Poor state, parameters need change!')
    else
        disp('improved success')
    end
end
M_noise = zeros(Nraw,Ncol);
M_signal = zeros(Nraw,Ncol);
for ii = 1:5
    if ifSignal(ii)
        M_signal = M_signal+u(:,:,ii);
    else
        M_noise = M_noise+u(:,:,ii);
    end
end
end

function [u, u_hat, omega,uDiff_set,omegaDiff_set] = VMD_2D(signal, alpha, tau, K, DC, init, tol,varargin)
% 2D Variational Mode Decomposition
% Authors: Konstantin Dragomiretskiy and Dominique Zosso
% {konstantin,zosso}@math.ucla.edu
% http://www.math.ucla.edu/~{konstantin,zosso}
% Initial release 2014-03-17 (c) 2014
%
% Input and Parameters:
% ---------------------
% signal     - the space domain signal (2D) to be decomposed
% alpha      - the balancing parameter for data fidelity constraint
% tau        - time-step of dual ascent ( pick 0 for noise-slack )
% K          - the number of modes to be recovered
% DC         - true, if the first mode is put and kept at DC (0-freq)
% init       - 0 = all omegas start at 0
%              1 = all omegas start initialized randomly
% tol        - tolerance of convergence criterion; typically around 1e-7
%
% When using this code, please do cite our papers:
% -----------------------------------------------
% K. Dragomiretskiy, D. Zosso, Variational Mode Decomposition, IEEE Trans.
% on Signal Processing, 62(3):531-544, 2014. DOI:10.1109/TSP.2013.2288675
%
% K. Dragomiretskiy, D. Zosso, Two-Dimensional Variational Mode
% Decomposition, IEEE Int. Conf. Image Proc. (submitted). Preprint
% available here: ftp://ftp.math.ucla.edu/pub/camreport/cam14-16.pdf
%


% Resolution of image
[Hy,Hx] = size(signal);
[X,Y] = meshgrid((1:Hx)/Hx, (1:Hy)/Hy);


% Spectral Domain discretization
fx = 1/Hx;
fy = 1/Hy;
freqs_1 = X - 0.5 - fx;
freqs_2 = Y - 0.5 - fy;

% N is the maximum number of iterations
N=300;

% For future generalizations: alpha might be individual for each mode
Alpha = alpha*ones(K,1);

% Construct f and f_hat
f_hat = fftshift(fft2(signal));

% Storage matrices for (Fourier) modes. All iterations are not recorded.
u_hat = zeros(Hy,Hx,K);
u_hat_old = u_hat;
sum_uk = 0;

% Storage matrices for (Fourier) Lagrange multiplier.
mu_hat = zeros(Hy,Hx);

% N iterations at most, 2 spatial coordinates, K clusters
omega = zeros(N, 2, K);

% Initialization of omega_k
switch init
    case 0
        % spread omegas radially
        
        % if DC, keep first mode at 0,0
        if DC
            maxK = K-1;
        else
            maxK = K;
        end
        for k = DC + (1:maxK) 
            omega(1,1,k) = 0.25*cos(pi*(k-1)/maxK);
            omega(1,2,k) = 0.25*sin(pi*(k-1)/maxK);
        end
        
        % Case 1: random on half-plane
    case 1
        for k=1:K
            omega(1,1,k) = rand()-1/2;
            omega(1,2,k) = rand()/2;
        end
        
        % DC component (if expected)
        if DC == 1
            omega(1,1,1) = 0;
            omega(1,2,1) = 0;
        end
    case 2
        omega(1,1,1:K) = varargin{1}(1,:);
        omega(1,2,1:K) = varargin{1}(2,:);
        % DC component (if expected)
        if DC == 1
            omega(1,1,1) = 0;
            omega(1,2,1) = 0;
        end
end

%% Main loop for iterative updates

% Stopping criteria tolerances
uDiff=tol+eps;
omegaDiff = tol+eps;
uDiff_set = zeros(N,1);
omegaDiff_set = zeros(N,1);
% first run
n = 1;

% run until convergence or max number of iterations
while( ( uDiff > tol && omegaDiff > tol ) && n < N )
% while( ( uDiff > tol || omegaDiff > tol ) && n < N )
    
    % first things first
    k = 1;
    
    % compute the halfplane mask for the 2D "analytic signal"
    HilbertMask = (sign(freqs_1*omega(n,1,k) + freqs_2*omega(n,2,k))+1);
    
    % update first mode accumulator
    sum_uk = u_hat(:,:,end) + sum_uk - u_hat(:,:,k);
    
    % update first mode's spectrum through wiener filter (on half plane)
    u_hat(:,:,k) = ((f_hat - sum_uk - mu_hat(:,:)/2).*HilbertMask)./(1+Alpha(k)*((freqs_1 - omega(n,1,k)).^2+(freqs_2 - omega(n,2,k)).^2));
    
    % update first mode's central frequency as spectral center of gravity
    if ~DC
        omega(n+1,1,k) = sum(sum(freqs_1.*(abs(u_hat(:,:,k)).^2)))/sum(sum(abs(u_hat(:,:,k)).^2));
        omega(n+1,2,k) = sum(sum(freqs_2.*(abs(u_hat(:,:,k)).^2)))/sum(sum(abs(u_hat(:,:,k)).^2));
        
        % keep omegas on same halfplane
        if omega(n+1,2,k) < 0
            omega(n+1,:,k) = -omega(n+1,:,k);
        end
    end
    
    % recover full spectrum from analytic signal
    u_hat(:,:,k) = fftshift(fft2(real(ifft2(ifftshift(squeeze(u_hat(:,:,k)))))));
    
    % work on other modes
    for k=2:K
        
        % recompute Hilbert mask
        HilbertMask = (sign(freqs_1*omega(n,1,k) + freqs_2*omega(n,2,k))+1);
        
        % update accumulator
        sum_uk = u_hat(:,:,k-1) + sum_uk - u_hat(:,:,k);
        
        % update signal spectrum
        u_hat(:,:,k) = ((f_hat - sum_uk - mu_hat(:,:)/2).*HilbertMask)./(1+Alpha(k)*((freqs_1 - omega(n,1,k)).^2+(freqs_2 - omega(n,2,k)).^2));
        
        % update signal frequencies
        omega(n+1,1,k) = sum(sum(freqs_1.*(abs(u_hat(:,:,k)).^2)))/sum(sum(abs(u_hat(:,:,k)).^2));
        omega(n+1,2,k) = sum(sum(freqs_2.*(abs(u_hat(:,:,k)).^2)))/sum(sum(abs(u_hat(:,:,k)).^2));
        
        % keep omegas on same halfplane
        if omega(n+1,2,k) < 0
            omega(n+1,:,k) = -omega(n+1,:,k);
        end
        
        % recover full spectrum from analytic signal
        u_hat(:,:,k) = fftshift(fft2(real(ifft2(ifftshift(squeeze(u_hat(:,:,k)))))));
    end
    
    % Gradient ascent for augmented Lagrangian
    mu_hat(:,:) = mu_hat(:,:) + tau*(sum(u_hat,3) - f_hat);
    
    % increment iteration counter
    n = n+1;
    
    % convergence?
    uDiff = eps;
    omegaDiff = eps;
    
    for k=1:K
        omegaDiff = omegaDiff + sum(sum(abs(omega(n,:,:) - omega(n-1,:,:)).^2));
        uDiff = uDiff + sum(sum(1/(Hx*Hy)*(u_hat(:,:,k)-u_hat_old(:,:,k)).*conj((u_hat(:,:,k)-u_hat_old(:,:,k)))));
    end
    
    uDiff = abs(uDiff);
    uDiff_set(n,1) = uDiff;
    omegaDiff_set(n,1) = omegaDiff;
    
    u_hat_old = u_hat;

end
fprintf('run %d times\n',n)

%% Signal Reconstruction

% Inverse Fourier Transform to compute (spatial) modes
u = zeros(Hy,Hx,K);
for k=1:K
    u(:,:,k) = real(ifft2(ifftshift(squeeze(u_hat(:,:,k)))));
end

% Should the omega-history be returned, or just the final results?
%omega = omega(n,:,:);
omega = omega(1:n,:,:);

end

function [xgd,ifSignal] = FFTA4VMD(u,varargin)
% resolution 3deg
%%
fft2d1 = fftshift(fft2(u(:,:,1)));
fft2d2 = fftshift(fft2(u(:,:,2)));
fft2d3 = fftshift(fft2(u(:,:,3)));
fft2d4 = fftshift(fft2(u(:,:,4)));
fft2d5 = fftshift(fft2(u(:,:,5)));
xgd = zeros(5,1);
ifSignal = nan(5,1);
% x = 1:61;
% y = 1:121;
% [ygrd,xgrd] = meshgrid(y,x);
if ~isempty(varargin)
    resolution = varargin{1};
    Nlon = 360/resolution;
else
    Nlon = size(u,2)-1;
    resolution = 360/Nlon;
end



x = 1:(Nlon+1);
y = gaussmf(x,[45/resolution/3 (1*Nlon/2)]) ;
lfilt = y(1:Nlon+1);
%%
% 取上界
% 利用mean+1std
% 对每一经度的求mean 和std\
% figure,
for ii = 1:5
%     eval(['input=abs(fft2d',num2str(ii),');']);
    switch ii
        case 1
            input=abs(fft2d1);
        case 2
            input=abs(fft2d2);
        case 3
            input=abs(fft2d3);
        case 4
            input=abs(fft2d4);
        case 5
            input=abs(fft2d5);
    end

lonm = zeros(Nlon+1,1);
lons = zeros(Nlon+1,1);
for i = 1:Nlon+1
    lonm(i) = mean(input(:,i));
    lons(i) = std(input(:,i));
end

lonupper = lonm+1*lons;
norm_lup = lonupper/(max(lonupper)-min(lonupper));
norm_lft = lfilt/(max(lfilt)-min(lfilt));
% xgd = corr(lonupper,lfilt')
% xgd = corr(norm_lup,norm_lft')
xgd(ii) = sum(norm_lup.*norm_lft');
% new criterian
[~,locs] = findpeaks(xcorr(norm_lft,norm_lup),1,'MinPeakProminence',0.5,'sortstr','descend');
[r,lags] = xcorr(norm_lft,norm_lup);
% ntitle([num2str(locs)])
if any(lags(locs)>=0-5&lags(locs)<=0+5)
    ifSignal(ii) = 1;
else
    ifSignal(ii) = 0;
end

end

end

