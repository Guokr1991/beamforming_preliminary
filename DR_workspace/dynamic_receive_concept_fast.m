function [focus_rf, x, z] = dynamic_receive_concept_fast(rf,acq_params,bf_params)
%testing for github

% rf - raw rf data organized [rf_line,rx_chan, tx_event]

%all focusing performed on axis

% acq_params.rx_pos; %lateral positions of each rx element relative to tx
% bf_params.x; %lateral positions of each tx focus (per tx event)
% acq_params.c;
% acq_params.t0; start time of rf data reference in terms of sample number
tic

x = bf_params.x;
% z_ref_test = single(acq_params.t0+1:size(rf,1)/acq_params.fs)*acq_params.c;
z_ref = ((acq_params.t0+1:acq_params.t0+size(rf,1))/acq_params.fs)*acq_params.c;
z = z_ref/2;

[b,a]=butter(2,[.05 .95]);

% intialize matrices and arrays for speed
n_tx = length(bf_params.x);
n_rcv_chn = length(acq_params.rx_pos);
n_depth = length(z);
focus_rf = zeros(length(z),n_tx);

if n_rcv_chn ~= size(rf,2) || n_tx ~= size(rf,3)
    disp('Mismatch in RF data.')
end
    
dz = repmat(z',1,n_rcv_chn);
dx = repmat(acq_params.rx_pos,n_depth,1);
dr = sqrt(dz.^2+dx.^2);
t_samp = (dr+dz)./acq_params.c;
    
%     ttrange=minmax(dz(:)');
%     trrange=minmax(dr(:)');
%     idx_l=find(z_ref<ttrange(1)+trrange(1),1,'last'); if(isempty(idx_l)), idx_l=1; end
%     idx_h=find(z_ref>ttrange(2)+trrange(2),1,'first'); if(isempty(idx_h)), idx_h=length(z_ref); end
    

%     focus_rf(:,j) = sum(lin_interp(z_ref(idx_l:idx_h)'/acq_params.c,squeeze(rf(idx_l:idx_h,:,j)),t_samp),2);
% rf_reshape = reshape(rf,[size(rf,1), size(rf,2)*size(rf,3)]);
% t_samp_rep = repmat(t_samp,1,n_tx);
% focus_rf = lin_interp(z_ref'/acq_params.c,rf_reshape,t_samp_rep);

for j = 1:n_tx % iterate for every tx event
    focus_rf(:,j) = sum(lin_interp(z_ref'/acq_params.c,squeeze(rf(:,:,j)),t_samp),2);
end

focus_rf = filter(b,a,focus_rf);
focus_rf(find(isnan(focus_rf))) = 0; 
toc

