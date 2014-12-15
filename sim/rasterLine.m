function [raster] = rasterLine(raster, x1, x2, y1, y2, v)
% Fill a line from (x1,y1) to (x2,y2) in a raster plot with the value v.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab seems to make a copy of the raster when calling to            %
% rasterLine: so it is much better to embed this function code in the  %
% calling code.                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dx = x2 - x1;
dy = y2 - y1;

if(dx || dy)
    if(abs(dx) >= abs(dy))
        if(dx>0)
            raster(sub2ind(size(raster), floor(linspace(y1, y2, dx+1)), x1:x2)) = v;
        else
            raster(sub2ind(size(raster), floor(linspace(y1, y2, 1-dx)), x1:-1:x2)) = v;
        end
    else
        if(dy>0)
            raster(sub2ind(size(raster), y1:y2, floor(linspace(x1, x2, dy+1)))) = v;
        else
            raster(sub2ind(size(raster), y1:-1:y2, floor(linspace(x1, x2, 1-dy)))) = v;
        end
    end
else
    raster(y1, x1) = v;
end
        
            