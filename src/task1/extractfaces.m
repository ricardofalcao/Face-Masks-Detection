function [BBx, New_BW] = extractfaces(BW)
    BBx = zeros([size(BW), 1]);
    I = size(BBx, 3);
    
    New_BW = BW;
    
    [L, N] = bwlabel(BW);
    
    for k = 1 : N
        Reg = (L == k);

        props = regionprops(Reg, 'BoundingBox', 'FilledImage', 'ConvexHull');
        filled = regionprops(props.FilledImage, 'Eccentricity', 'MajorAxisLength', 'MinorAxisLength');

        if (0.7 < filled.Eccentricity && filled.Eccentricity < 0.9) && (filled.MajorAxisLength < 2.5 * filled.MinorAxisLength)
            fprintf('Object has Eccentricity = %f, AxisRatio = %f\n', filled.Eccentricity, filled.MajorAxisLength / filled.MinorAxisLength);

            BB = props.BoundingBox;
            BB_BW = zeros(size(L));

            for x = 0 : uint16(BB(3))
                for y = 0 : uint16(BB(4))

                    BB_BW(uint16(BB(2)) + y, uint16(BB(1)) + x) = 1;
                end
            end

            New_BW = New_BW & ~(Reg & BB_BW);  
           
            BBx(:,:,I) = BB_BW;
            I = I + 1;
        end

    end

end