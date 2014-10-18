close all; clear all; clc;

addpath('/home/wjl11/Documents/MATLAB/field_ii');
addpath('/home/wjl11/Documents/MATLAB/field_ii_supplementary');

PARAM.DR = 1;
    PARAM.DR_multi = 0;
PARAM.plane = 0;
    PARAM.saveRF = 0;

% define transducer
f0 = 5e6;
BW = 0.7;
N_el = 128;
elv_focus = 0.04;
elv_fnum = 6.5;
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

% define phantom at different scanning angles
r = 3; % scan radius in cm

pr_on = [-1/3 0 1/3];
pr_off = [repmat([-0.5 1]',[3 1]),...
    [repmat(2,[2 1]); repmat(3, [2 1]); repmat(4, [2 1])]];

x = pr_off(:,1);
z = r-pr_off(:,2);
th = [-20:20:20].*pi/180; % 0 at center, scanning perpendicular to phantom
for i = 1:length(th)
    th_i = atan(x./z);
    h = x./sin(th_i);
    phx_off(i,:) = -h.*sin(th(i)-th_i);
    phz_off(i,:) = r-h.*cos(th(i)-th_i);
end

for i = 1:length(th)
    phx_on(i,:) = -pr_on*r*sin(th(i));
    phz_on(i,:) = r*(1-pr_on*cos(th(i)));
end

ph_x = [phx_off phx_on];
ph_z = [phz_off phz_on];

figure; hold on
for i = 1:(length(pr_on)+length(pr_off))
    plot(ph_x(:,i),ph_z(:,i),'x','Linewidth',2)
end
hold off

% intialize point target configurations
for i = 1:length(th)
    position_pts = [ph_x(i,:)',zeros(length(pr_on)+length(pr_off),1), ph_z(i,:)'].*1e-2;
    amplitude_pts = 100.*ones(length(pr_on)+length(pr_off),1);
    
    if PARAM.DR == 1
        x = -0.02:0.0005:0.02;
        foc_z = 0.03; 
        if PARAM.DR_multi == 1 % walk aperture for calc_scat_multi DR
            for nn = 1:length(x) 
                xdc_center_focus(tx,[x(nn) 0 0]);
                xdc_focus(tx,0,[x(nn) 0 foc_z]);
                xdc_center_focus(rx,[x(nn) 0 0]);
                xdc_dynamic_focus(rx,0,0,0);
                [tmp_rf{nn}, tmp_st(nn)] = calc_scat_multi(tx, rx, position_pts, amplitude_pts);
            end
            [pad_rf s0 t0] = shift_times_multi(tmp_rf, tmp_st,fs);
            rf = im_sum(pad_rf,3);
        elseif PARAM.DR_multi == 0 % walk aperture for calc_scat DR
            for nn = 1:length(x) 
                xdc_center_focus(tx,[x(nn) 0 0]);
                xdc_focus(tx,0,[x(nn) 0 foc_z]);
                xdc_center_focus(rx,[x(nn) 0 0]);
                xdc_dynamic_focus(rx,0,0,0);
                [temp, start(nn)] = calc_scat(tx,rx,position_pts, amplitude_pts);
                rf(1:length(temp), nn) = temp;
            end
            [rf, t0]=ShiftStartTimes(start, rf, fs);
        end
        
        z = (t0+(1:size(rf,1))/fs)*1540/2;
        env = abs(hilbert(rf));
        env_db = 20*log10(env/max(env(:)));
        figure;
        imagesc(x,z,env_db,[-60 0]); axis image; colormap gray;
        xlabel('x (cm)'),ylabel('z (cm)');
        title(['Image at ' num2str(th(i)*180/pi) ' degrees'])
        clear rf start
    end
    
    if PARAM.plane == 1 % plane wave tx for calc_scat_multi SSA
        Ncyc = 2;
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
        v = squeeze(rf(:,:,i));
        v = v./max(max(v));
        t = (0:size(rf,1)-1)/fs+t0;
%         figure; hold on
%         for j = 1:N_el
%             plot(v(:,j)+j, t);
%         end
%         hold off
%         set(gca,'YDir','reverse')
%         title(['RF at ' num2str(th(i)*180/pi) ' degrees'])
%         grid on
%         xlim([1 128])
    end
    
    if PARAM.saveRF == 1
        acq_params.r = r;
        acq_params.th = th;
        acq_params.c = c;
        acq_params.fs = fs;
        acq_params.f0 = f0;
        acq_params.s0 = s0;
        acq_params.rx_pos = rx_pos;
        save('/home/wjl11/Documents/MATLAB/beamforming_preliminary/SSA_workspace/field_sim/SSA_field_test2.mat','rf','acq_params');
    end
end

xdc_free(rx)
xdc_free(tx)
