function i = inside(d, w, X) % .................................... inside
% INSIDE - return true (false) for points inside (outside) a domain
%
% i = INSIDE(d, p) returns a logical array of values specifying if
%   each point in the list of complex numbers p is inside (1) or outside
%   (0) the domain d. The output i is the same shape as input p.
%
% Issues/notes:
%  * Uses fast inpolygon implementations in inpolywrapper
%  * inluded temporary hack to make semi-infinite strip domain if corners
%    don't connect up - the domain is officially borken then anyway!
%         if d.exterior
%           i = logical(ones(size(X(:))));
%         else                                % interior domain
%           js = find(d.spiece==0);           % indices of segs outer bdry piece
%           v = domain.approxpolygon(d.seg(js), d.pm(js));
v = w.';
i = inpolywrapper(X(:), v);
% end
%         for piece=1:max(d.spiece)         % kill pts from each interior piece
%           js = find(d.spiece==piece);
%           v = domain.approxpolygon(d.seg(js), d.pm(js));
%           % HACKs which extend an excluded strip down/upwards: (all pm equal)
%           % (this is used in quasi-periodic scattering situations)
%           if isnan(d.cloc)
%             if d.pm(js(1))==-1  % usual way around; exclude domain below
%               x0 = d.seg(js(end)).eloc(1); x1 = d.seg(js(1)).eloc(2);
%               v = [x1-10i; v; x0; x0-10i]; % works for multiple segs
%             else                % exclude domain above (original code)
% %              e=(d.pm(js)+3)/2; v = [d.seg(js).eloc(3-e)+d.pm(js)*10i; v; d.seg(js).eloc(e); d.seg(js).eloc(e)+d.pm(js)*10i];
%               x0 = d.seg(js(1)).eloc(1); x1 = d.seg(js(end)).eloc(2);
%               v = [x0+10i; v; x1; x1+10i];  % new code 3/16/13, multiple segs
% %              figure; plot([v;v(1)], 'k--'); % debug excluded region
%             end
%           end
%           i = i & ~utils.inpolywrapper(X(:), v);
%         end
i = reshape(i, size(X));
end