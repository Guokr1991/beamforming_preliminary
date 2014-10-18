% planewave.m - Steered plane wave transmits from center of the array
%
% check:
% acq_params
%   rx_pos - Receive channel positions, x or [x y z] one row per element
%            in meters
%   theta - Angles of each transmit event (degrees)
% bf_params
%   x,z - Grid to beamform final image on (m)
%   tx_channel - Preserve transmit channel data instead of receive
%   mask_radius - Perpendicular radius from center of plane wave to
%                 preserve for each transmit event
function varargout = planewave(mode,acq_params,bf_params,rf)

    switch mode
        case 'check_params'
            [acq_params,bf_params]=check_params(acq_params,bf_params);
            varargout{1}=acq_params;
            varargout{2}=bf_params;
        case 'calc_delays'
            [tt tr mask]=calc_delays(acq_params,bf_params);
            varargout{1}=tt;
            varargout{2}=tr;
            varargout{3}=mask;
        case 'pre_focus'
            varargout{1}=pre_focus(acq_params,bf_params,rf);
        case 'post_focus'
            varargout{1}=post_focus(acq_params,bf_params,rf);
    end

end

function [acq_params,bf_params] = check_params(acq_params,bf_params)
    if(size(acq_params.rx_pos,2)~=3), acq_params.rx_pos=[acq_params.rx_pos(:),zeros(length(acq_params.rx_pos),2)]; end
    if(~isfield(bf_params,'z')), bf_params.z=single((acq_params.t0+1:acq_params.t0+acq_params.samples)/acq_params.fs)*acq_params.c/2; bf_params.z=bf_params.z(bf_params.z>0); end
    if(~isfield(bf_params,'x')), bf_params.x=single(acq_params.rx_pos(:,1)); end
    bf_params.y=single(0);
    if(~isfield(bf_params,'tx_channel')), bf_params.tx_channel=0; end
    if(~isfield(bf_params,'mask_radius')), bf_params.mask_radius=Inf; end
end

function [tt tr mask] = calc_delays(acq_params,bf_params)
    [z,x]=ndgrid(bf_params.z,bf_params.x);
    x=x(:);z=z(:);
       
    tt=zeros(length(x),acq_params.txEvents,'single');
    for i=1:acq_params.txEvents
        cur_focus=[cosd(-acq_params.theta(i)) sind(-acq_params.theta(i))];
        tt(:,i)=abs(cur_focus(2)*x-cur_focus(1)*z)/sqrt(sum(cur_focus.^2));
    end
    
    tr=zeros(length(x),acq_params.rxChannels,'single');
    for i=1:size(acq_params.rx_pos,1)
        tr(:,i)=sqrt((acq_params.rx_pos(i,1)-x).^2+(acq_params.rx_pos(i,2)).^2+(acq_params.rx_pos(i,3)-z).^2);
    end
    
    if(bf_params.tx_channel)     
        %Swap transmit/receive roles
        temp=tr;
        tr=tt;
        tt=temp;
        mask=true(length(x),acq_params.rxChannels);
    else
        mask=false(length(x),acq_params.txEvents);
        for i=1:acq_params.txEvents
            cur_focus=[sind(acq_params.theta(i)) cosd(acq_params.theta(i))];
            mask(:,i)=abs(cur_focus(2)*x-cur_focus(1)*z)/sqrt(sum(cur_focus.^2))<bf_params.mask_radius;
        end
    end
    
    
end

function rf = pre_focus(acq_params,bf_params,rf_in)
    if(bf_params.tx_channel)
        rf=permute(rf_in,[1 3 2]);
    else
        rf=rf_in;
    end
end

function rf = post_focus(acq_params,bf_params,rf_in)
    rf=rf_in;
    if(bf_params.tx_channel)
        [z,x]=ndgrid(bf_params.z,bf_params.x);
        x=x(:);z=z(:);
        for i=1:acq_params.txEvents
            cur_focus=[sind(acq_params.theta(i)) cosd(acq_params.theta(i))];
            mask=abs(cur_focus(2)*x-cur_focus(1)*z)/sqrt(sum(cur_focus.^2))>=bf_params.mask_radius;
            rf(mask,i)=0;
        end
    end
end