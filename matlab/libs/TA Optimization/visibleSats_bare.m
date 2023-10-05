function [visible] = visibleSats_bare(r, sats, time)
%VISIBLESATS_BARE Compute number of satellites visible to the user (on
%lunar surface).
%   Stripped down version of visibleSats(); only considers a user at the
%   south pole with 0 deg elevation angle.

visible = [];
for i=1:size(sats,3)
    if sats(time,3,i) < -r
        visible = [visible i];
    end
end
end

