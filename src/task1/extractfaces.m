function [BBx, NBB, New_BW, L] = extractfaces(BW)
    Debug = 1;
    
    Size = size(BW);
    BBx = zeros([Size, 1]);
    NBB = size(BBx, 3) - 1;
    
    New_BW = BW;
    
    [L, N] = bwlabel(BW);
    
    for k = 1 : N
        Reg = (L == k);
                
        if Debug == 1
            fprintf('Region %d\n', k);
        end

        props = regionprops(Reg, 'BoundingBox', 'FilledImage', 'ConvexImage', 'Area', 'Solidity', 'Centroid');
        BB = props.BoundingBox;
        C = props.Centroid;
        
        NullThickness = 0.05;
        if (C(1) < NullThickness*Size(1) || C(1) > (1 - NullThickness)*Size(1))
            continue
        end
        
        if (C(2) < NullThickness*Size(2) || C(2) > (1 - NullThickness)*Size(2))
            continue
        end
        
        if (props.Solidity < 0.7)
            if Debug == 1
                fprintf('Object failed Solidity -> %f\n', props.Solidity);
            end
            
            continue;
        end       
        
        
        if (props.Area < 3000)
            if Debug == 1
                fprintf('Object failed Area -> %f\n', props.Area);
            end
            
            New_BW = imabsdiff(New_BW, Reg);
            continue;
        end

        WHRatio = BB(3) / BB(4);
        
        if (WHRatio < 0.4 || WHRatio > 1.8)
            if Debug == 1
                fprintf('Object failed Width to Height Ratio -> %f\n', WHRatio);
            end
            
            continue;
        end
        
        filled = regionprops(props.FilledImage, 'Eccentricity', 'Orientation', 'MajorAxisLength', 'MinorAxisLength');
        convex = regionprops(props.ConvexImage, 'Eccentricity');
        
        o = abs(filled.Orientation);
        if (o <= 50 || o >= 130)
            if Debug == 1
                fprintf('Object failed Orientation -> %f\n', filled.Orientation);
            end
            
            continue;
        end
        
        if (filled.Eccentricity <= 0.7 || filled.Eccentricity >= 0.9)
            if Debug == 1
                fprintf('Object failed Eccentricity -> %f\n', filled.Eccentricity);
            end
            
            continue;
        end
        
        if (abs(filled.Eccentricity - convex.Eccentricity) > 0.07)
            if Debug == 1
                fprintf('Object failed absolute Eccentricity -> %f\n', abs(filled.Eccentricity - convex.Eccentricity));
            end
            
            continue;
        end
        
        if Debug == 1
            fprintf('Object has Orientation = %f, Eccentricity = %f, Convex.Eccentricity = %f, WHRatio = %f, Area = %f\n', filled.Orientation, filled.Eccentricity, convex.Eccentricity, WHRatio, props.Area);
        end
        
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