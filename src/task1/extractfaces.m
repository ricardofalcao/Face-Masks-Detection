function [BBx, NBB, New_BW, new_maxArea] = extractfaces(BW, maxArea)
    Debug = 0;
    TP_info = 0;
    
    new_maxArea = maxArea;
    
    Size = size(BW);
    BBx = zeros([Size, 1]);
    NBB = size(BBx, 3) - 1;
    
    New_BW = BW;
    
    [L, N] = bwlabel(BW);
    
    for k = 1 : N
        Reg = (L == k);

        props = regionprops(Reg, 'BoundingBox', 'FilledImage', 'ConvexImage', 'Area', 'Solidity', 'Circularity', 'Centroid');
        BB = props.BoundingBox;
        C = props.Centroid;
        
        Round = props.Circularity > 0.7;
                
        if Debug == 1
            fprintf('Region %d:\n', k, props.Circularity);
        end
        
        NullThickness = 0.01;
        if (C(1) < NullThickness*Size(1) || C(1) > (1 - NullThickness)*Size(1))
            if Debug == 1
                fprintf('Object failed: Outside valid region\n');
            end
            
            continue
        end
        
        if (C(2) < NullThickness*Size(2) || C(2) > (1 - NullThickness)*Size(2))
            if Debug == 1
                fprintf('Object failed: Outside valid region\n');
            end
            
            continue
        end
        
        if (props.Solidity < 0.6 || props.Solidity > 0.96)
            if Debug == 1
                fprintf('Object failed Solidity -> %f\n', props.Solidity);
            end
            
            continue;
        end       
        
         Rectangularity = props.Area / (BB(3)*BB(4));
         
         if (~Round && (Rectangularity < 0.5 || Rectangularity > 0.77))
            if Debug == 1
                fprintf('Object failed Rectangularity -> %f\n', Rectangularity);
            end
            
            continue;
        end          
        
        if (props.Area < 5100 || props.Area > 365000)
            if Debug == 1
                fprintf('Object failed Area -> %f\n', props.Area);
            end
            
            if (props.Area < 5000)
                New_BW = imabsdiff(New_BW, Reg);
            end
            
            continue;
        end

        WHRatio = BB(3) / BB(4);
        
        if (WHRatio < 0.45 || WHRatio > 1.2)
            if Debug == 1
                fprintf('Object failed Width to Height Ratio -> %f\n', WHRatio);
            end
            
            continue;
        end
        
        filled = regionprops(props.FilledImage, 'Eccentricity', 'MajorAxisLength', 'MinorAxisLength');
        convex = regionprops(props.ConvexImage, 'Eccentricity', 'Orientation');
        
        o = abs(convex.Orientation);
        if (~Round && o <= 15)
            if Debug == 1
                fprintf('Object failed Orientation -> %f\n', convex.Orientation);
            end
            
            continue;
        end

        if (~Round && (filled.Eccentricity <= 0.6 || filled.Eccentricity >= 0.92))
            if Debug == 1
                fprintf('Object failed Eccentricity -> %f\n', filled.Eccentricity);
            end
            
            continue;
        end
        
        if (~Round && (abs(filled.Eccentricity - convex.Eccentricity) > 0.1))
            if Debug == 1
                fprintf('Object failed absolute Eccentricity -> %f\n', abs(filled.Eccentricity - convex.Eccentricity));
            end
            
            continue;
        end
        
        if maxArea > 4 * props.Area
             if Debug == 1
                fprintf('Object failed Relative Area -> %f > 4* %f\n', maxArea, props.Area);
             end
             
            New_BW = imabsdiff(New_BW, Reg);
            
            continue;
           
        end
        
        if TP_info == 1
            fprintf('Object has Orientation = %f, Eccentricity = %f, Convex.Eccentricity = %f, WHRatio = %f, Area = %f\n', convex.Orientation, filled.Eccentricity, convex.Eccentricity, WHRatio, props.Area);
            fprintf('Object has Solidity = %f, Rectangularity = %f\n', props.Solidity, Rectangularity);
        end
        
        BB_BW = zeros(size(L));

        for x = 1 : floor(BB(3))
            for y = 1 : floor(BB(4))

                BB_BW(floor(BB(2)) + y, floor(BB(1)) + x) = 1;
            end
        end
        
        New_BW = New_BW & ~(Reg & BB_BW);  

        BBx(:,:,NBB + 1) = BB_BW;
        NBB = NBB + 1;
        
        new_maxArea = max([nnz(Reg), maxArea]);

        
    end

end