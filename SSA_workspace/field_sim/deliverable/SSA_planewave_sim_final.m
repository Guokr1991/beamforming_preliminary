close all; clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% RUN CONIGURATIONS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PARAM.DR = 0;
    PARAM.DR_multi = 0;
PARAM.plotPh = 1;
PARAM.plane = 1;
    PARAM.saveRF = 1;
    PARAM.plotRF = 1;
    
% define phantom at different scanning angles
scan_r = 6;
gridx = [-1 0 1 0]; % x point coord relative to focal center [cm]
gridz = [1 1 -1 0];
ph_pts(:,:,1) = repmat(gridx,[length(gridz) 1]);
ph_pts(:,:,2) = repmat(gridz',[1 length(gridx)]);
th = (-40:5:40).*pi/180;
if length(th)>4
    PARAM.DR=0;
    PARAM.plotRF=0;
    disp('Number of frames exceed limit - no plot generated.')
end
    
% add field paths and field_init    
addpath('/home/wjl11/Documents/MATLAB/field_ii');
addpath('/home/wjl11/Documents/MATLAB/field_ii_supplementary');
field_init(-1);

% define transducer
f0 = 5e6;
BW = 0.7;
N_el = 128;
elv_focus = 2; % CHANGED FOR PLANE WAVE w/o ELECTRONIC FOCUSING
elv_fnum = 6.5*elv_focus/0.04;
kerf_fraction = 0.05;
focus = [0 0 elv_focus];

% general parameters
c = 1540; fs = 100e6;

% derived parameters
lambda = c/f0;
el_height = elv_focus/elv_fnum;
el_width = (1-kerf_fraction)*lambda;
el_kerf = kerf_fraction*lambda;
n_sub_x = ceil(el_width/(lambda/2));
n_sub_y = ceil(el_height/(lambda/2));
el_pitch=el_width+el_kerf;
rx_pos = (-((N_el/2-1)*el_pitch+el_pitch/2):el_pitch:((N_el/2-1)*el_pitch+el_pitch/2));

% create transmit and receive arrays
tx=xdc_focused_array(N_el,el_width,el_height,el_kerf,elv_focus,n_sub_x,n_sub_y,focus);
rx=xdc_focused_array(N_el,el_width,el_height,el_kerf,elv_focus,n_sub_x,n_sub_y,focus);

% set impulse responses
tc=gauspuls('cutoff',f0,BW,-6,-40);
imp_resp=gauspuls(-tc:1/fs:tc,f0,BW);
xdc_impulse(tx,imp_resp);
xdc_impulse(rx,imp_resp);

% define points in phantom
r = scan_r/100;
x = ph_pts(:,:,1)./100;
z = r - ph_pts(:,:,2)./100;

xi_tmp = zeros(size(x));
zi_tmp = zeros(size(z));
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

    xi_tmp = double(xi_tmp);
    zi_tmp = double(zi_tmp);

    ph_x(:,nf) = xi_tmp(:);
    ph_z(:,nf) = zi_tmp(:);
end

if PARAM.plotPh
    i = find(th==0);
    figure;
    plot(ph_x(:,i), ph_z(:,i),'o','Linewidth',2);
    set(gca,'YDir','reverse')
    grid on;
    title(['Focal points for frame at ' num2str(th(i)*180/pi) ' degrees'])
end

% intialize point target configurations
for i = 1:length(th)
    position_pts = [ph_x(:,i),zeros(length(ph_x(:,i)),1), ph_z(:,i)];
    amplitude_pts = 100.*ones(length(ph_x(:,i)),1);
 
    if PARAM.DR == 1
        x_DR = -0.02:0.0005:0.02; 
        foc_z = 0.03; 
        if PARAM.DR_multi == 1 % walk aperture for calc_scat_multi DR
            for nn = 1:length(x_DR) 
                xdc_center_focus(tx,[x_DR(nn) 0 0]);
                xdc_focus(tx,0,[x_DR(nn) 0 foc_z]);
                xdc_center_focus(rx,[x_DR(nn) 0 0]);
                xdc_dynamic_focus(rx,0,0,0);
                [tmp_rf{nn}, tmp_st(nn)] = calc_scat_multi(tx, rx, position_pts, amplitude_pts);
            end
            [pad_rf s0 t0] = shift_times_multi(tmp_rf, tmp_st,fs);
            rf = sum(pad_rf,3);
        elseif PARAM.DR_multi == 0 % walk aperture for calc_scat DR
            for nn = 1:length(x_DR) 
                xdc_center_focus(tx,[x_DR(nn) 0 0]);
                xdc_focus(tx,0,[x_DR(nn) 0 foc_z]);
                xdc_center_focus(rx,[x_DR(nn) 0 0]);
                xdc_dynamic_focus(rx,0,0,0);
                [temp, start(nn)] = calc_scat(tx,rx,position_pts, amplitude_pts);
                rf(1:length(temp), nn) = temp;
            end
            [rf, t0]=ShiftStartTimes(start, rf, fs);
        end
        
        z_DR = (t0+(1:size(rf,1))/fs)*1540/2;
        env = abs(hilbert(rf));
        env_db = 20*log10(env/max(env(:)));
        figure;
        imagesc(x_DR,z_DR,env_db,[-60 0]); axis image; colormap gray;
        xlabel('x (cm)'),ylabel('z (cm)');
        title(['Image at ' num2str(th(i)*180/pi) ' degrees'])
        clear rf start
    end
    
    if PARAM.plane == 1 % plane wave tx for calc_scat_multi SSA
        Ncyc = 1;
        t_ex = 0:1/fs:Ncyc/f0;
        pulse_ex = sin(2*pi*f0*t_ex);
        xdc_excitation(tx, pulse_ex);
        [tmp_rf{i}, tmp_st(i)] = calc_scat_multi(tx, rx, position_pts, amplitude_pts);
    else
    end
end

if PARAM.plane == 1
    [rf s0 t0] = shift_times_multi(tmp_rf, tmp_st, fs);
    for i = 1:length(th)

        if PARAM.plotRF == 1
            show_RF_rcv(squeeze(rf(:,:,i)),['RF at ' num2str(th(i)*180/pi) ' degrees'],fs,t0);
        end
    end
    
    if PARAM.saveRF == 1
        acq_params.r = scan_r;
        acq_params.th = th;
        acq_params.c = c;
        acq_params.fs = fs;
        acq_params.f0 = f0;
        acq_params.s0 = s0;
        acq_params.rx_pos = rx_pos;
        str = input('Input simulation file name: ','s');
        save(['/home/wjl11/Documents/MATLAB/beamforming_preliminary/SSA_workspace/field_sim/' str '.mat'],'rf','acq_params');
        disp([str '.mat saved.'])
    end
end

xdc_free(rx)
xdc_free(tx)
clear all;