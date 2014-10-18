clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUT PARAMETERS %%%
%%%%%%%%%%%%%%%%%%%%%%%%
c = 1540; 
rx_pos = -0.01:0.005:0.01; % center of all rcv elements
scan_r = 3; % scan radius in [cm]
th = (-20:20:20).*pi/180;
gridx = -10:10; % x focal coord relative to focal center [cm]
gridz = -3:4; % z focal coord relative to focal center [cm]
focal_pts(:,:,1) = repmat(gridx,[length(gridz) 1]);
focal_pts(:,:,2) = repmat(gridz',[1 length(gridx)]);


%%%%%%%%%%%%%%%%%%%%%%%%
%%%    MAIN CODE     %%%
%%%%%%%%%%%%%%%%%%%%%%%%

r = scan_r/100;
x = (focal_pts(:,:,1))./100;
z = (r + focal_pts(:,:,2))./100;

xi_tmp = zeros(size(x));
zi_tmp = zeros(size(z));
xi = zeros(size(x,1), size(x,2), length(th));
zi = zeros(size(z,1), size(z,2), length(th));

for nf = 1:length(th) 
    th_i = atan(x./(r-z));
    % off axis focal points
    xi_tmp(th_i ~= 0) = -(x(th_i ~= 0)./sin(th_i(th_i ~= 0))).*sin(th(nf)-th_i(th_i ~= 0));
    zi_tmp(th_i ~= 0) = r-(x(th_i ~= 0)./sin(th_i(th_i ~= 0))).*cos(th(nf)-th_i(th_i ~= 0));
    % on axis focal points
    xi_tmp(th_i == 0) = -(r-z(th_i == 0)).*sin(th(nf));
    zi_tmp(th_i == 0) = r - (r-z(th_i == 0)).*cos(th(nf));
    % at scan center (z = r and x = 0)
    xi_tmp(isnan(th_i)) = 0;
    zi_tmp(isnan(th_i)) = r;

    xi(:,:,nf) = single(xi_tmp);
    zi(:,:,nf) = single(zi_tmp);
    
%     figure
%     plot(xi_tmp(:), zi_tmp(:),'x','Linewidth',2);
%     set(gca,'YDir','reverse')
%     axis image
%     title(['Focal points for frame at ' num2str(th(nf)*180/pi) ' degrees'])

    for nr = 1:length(rx_pos)
        ri = sqrt((rx_pos(nr)-xi(:,:,nf)).^2+zi(:,:,nf).^2);
        t_samp = (ri+zi(:,:,nf))./c;
%         lin_interp()
    end
end



