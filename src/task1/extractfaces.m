function [BBx, NBB, New_BW] = extractfaces(BW)
    BBx = zeros([size(BW), 1]);
    NBB = size(BBx, 3) - 1;
    
    New_BW = BW;
    
    [L, N] = bwlabel(BW);
    
    for k = 1 : N
        Reg = (L == k);

        props = regionprops(Reg, 'BoundingBox', 'FilledImage', 'ConvexImage', 'Area');
        BB = props.BoundingBox;
        
        if (props.Area < 3000)
            continue;
        end

        WHRatio = BB(3) / BB(4);
        
        if (WHRatio < 0.4 || WHRatio > 1.0)
            continue;
        end
        
        filled = regionprops(props.FilledImage, 'Eccentricity', 'Orientation', 'MajorAxisLength', 'MinorAxisLength');
        convex = regionprops(props.ConvexImage, 'Eccentricity');
        
        o = abs(filled.Orientation);
        if (o <= 70 || o >= 110)
            continue;
        end
        
        if (filled.Eccentricity <= 0.7 || filled.Eccentricity >= 0.92)
            continue;
        end
        
        fprintf('Object has Orientation = %f, Eccentricity = %f; %f, WHRatio = %f, Area = %f\n', filled.Orientation, filled.Eccentricity, convex.Eccentricity, 0, props.Area);
        
        BB_BW = zeros(size(L));

        for x = 0 : uint16(BB(3))
            for y = 0 : uint16(BB(4))

                BB_BW(uint16(BB(2)) + y, uint16(BB(1)) + x) = 1;
            end
        end

        New_BW = New_BW & ~(Reg & BB_BW);  

        BBx(:,:,NBB + 1) = BB_BW;
        NBB = NBB + 1;

    end

end