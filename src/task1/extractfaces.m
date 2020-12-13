function [BBx, NBB, New_BW, L, new_maxArea] = extractfaces(BW, maxArea)
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
        
        if (props.Solidity < 0.65 || props.Solidity > 0.92)
            if Debug == 1
                fprintf('Object failed Solidity -> %f\n', props.Solidity);
            end
            
            continue;
        end       
        
         Rectangularity = props.Area / (BB(3)*BB(4));
         
         if (Rectangularity > 0.8)
            if Debug == 1
                fprintf('Object failed Rectangularity -> %f\n', Rectangularity);
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
        
        if (WHRatio < 0.45 || WHRatio > 1)
            if Debug == 1
                fprintf('Object failed Width to Height Ratio -> %f\n', WHRatio);
            end
            
            continue;
        end
        
        filled = regionprops(props.FilledImage, 'Eccentricity', 'Orientation', 'MajorAxisLength', 'MinorAxisLength');
        convex = regionprops(props.ConvexImage, 'Eccentricity');
        
        o = abs(filled.Orientation);
        if (o <= 40)
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
        
        fprintf('Object has Orientation = %f, Eccentricity = %f, Convex.Eccentricity = %f, WHRatio = %f, Area = %f\n', filled.Orientation, filled.Eccentricity, convex.Eccentricity, WHRatio, props.Area);
        fprintf('Object has Solidity = %f, Rectangularity = %f\n', props.Solidity, Rectangularity);
        
        BB_BW = zeros(size(L));

        for x = 0 : uint16(BB(3))
            for y = 0 : uint16(BB(4))

                BB_BW(uint16(BB(2)) + y, uint16(BB(1)) + x) = 1;
            end
        end

        New_BW = New_BW & ~(Reg & BB_BW);  

        BBx(:,:,NBB + 1) = BB_BW;
        NBB = NBB + 1;
        
        new_maxArea = nnz(Reg);

    end

end