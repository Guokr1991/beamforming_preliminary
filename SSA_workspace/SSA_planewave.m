
function [focus_rf, xdim, zdim] = SSA_planewave(rf, acq_params)
% acq_params
%     th - aperture angles [rad]

addpath('/home/wjl11/Documents/MATLAB/beamforming_preliminary/common/')
%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUT PARAMETERS %%%
%%%%%%%%%%%%%%%%%%%%%%%%

c = acq_params.c; 
rx_pos = acq_params.rx_pos; % center of all rcv elements
scan_r = acq_params.r; % scan radius in [cm]
r = scan_r/100;
th = acq_params.th;
fs = acq_params.fs;
f0 = acq_params.f0;

wl = c/f0;
b = 2*max(rx_pos);
ax_res = 100* (1/fs)*c/2; % axial resolution converted [cm]
lat_res = 100*(wl*r/b)/2; % lateral resolution converted [cm]

rangex = [-1.1 1.1]; % focal range of x [cm]
rangez = [-1.1 1.1]; % focal range of z relative to scan rad [cm]

if length(rangex) == 1 && length(rangez) == 1
    gridx = rangex;
    gridz = rangez;
else
    gridx = lat_res.*(ceil(rangex(1)/lat_res):1:floor(rangex(2)/lat_res));
    gridz = ax_res.*(ceil(rangez(1)/ax_res):1:floor(rangez(2)/ax_res));
end
% gridx = -1:0.005:0; % x focal coord relative to focal center [cm]
% gridz = -1.5:0.005:0; % z focal coord relative to focal center [cm]

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%    MAIN    %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%


t_ref = (1:acq_params.s0+size(rf,1))/acq_params.fs;
% t_ref = ((acq_params.s0+1:acq_params.s0+size(rf,1))/acq_params.fs);
focal_pts(:,:,1) = repmat(gridx,[length(gridz) 1]);
focal_pts(:,:,2) = repmat(gridz',[1 length(gridx)]);

x = focal_pts(:,:,1)./100;
z = r - focal_pts(:,:,2)./100;

xi_tmp = zeros(size(x));
zi_tmp = zeros(size(z));
xi = zeros(size(x,1), size(x,2), length(th));
zi = zeros(size(z,1), size(z,2), length(th));

focus_rf = zeros(length(gridz),length(gridx),length(th));
for nt = 1:length(th)  
    th_i = atan(x./(r-z));
    % off axis focal points
    xi_tmp(th_i ~= 0) = -(x(th_i ~= 0)./sin(th_i(th_i ~= 0))).*sin(th(nt)-th_i(th_i ~= 0));
    zi_tmp(th_i ~= 0) = r-(x(th_i ~= 0)./sin(th_i(th_i ~= 0))).*cos(th(nt)-th_i(th_i ~= 0));
    % on axis focal points
    xi_tmp(th_i == 0) = -(r-z(th_i == 0)).*sin(th(nt));
    zi_tmp(th_i == 0) = r - (r-z(th_i == 0)).*cos(th(nt));
    % at scan center (z = r and x = 0)
    xi_tmp(isnan(th_i)) = 0;
    zi_tmp(isnan(th_i)) = r;

    xi(:,:,nt) = double(xi_tmp);
    zi(:,:,nt) = double(zi_tmp);
    
    foc_tmp = zeros(length(gridz), length(gridx), length(rx_pos));
    
    for nr = 1:length(rx_pos)
        ri = sqrt((rx_pos(nr)-xi(:,:,nt)).^2+zi(:,:,nt).^2);
        t_samp = (ri+zi(:,:,nt))./c;
        foc_tmp(:,:,nr) = lin_interp(t_ref',[zeros(acq_params.s0,1); rf(:,nr, nt)], t_samp); 
    end
    foc_tmp(isnan(foc_tmp)) = 0;
    focus_rf(:,:,nt) = sum(foc_tmp,3);
    disp(['Frame at ' num2str(th(nt)*180/pi) ' degrees processed.'])
end
    
xdim = gridx;
zdim = r*100-gridz;


