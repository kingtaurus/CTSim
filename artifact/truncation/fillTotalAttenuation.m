function totalAttenuations = fillTotalAttenuation( truncatedView,  totalAttenuations )
% function totalAttenuations = fillTotalAttenuation( truncatedView,  totalAttenuations )
% Detecting truncation projection in sinogram
%
% Meng Wu at Stanford University
% 2014

noViews = length( truncatedView );

for iv = 1 : noViews
    
    if truncatedView( iv )
        
        leftTA = 0;
        rightTA = 0;
        
        il = iv;
        while il > 1 
            il = il - 1;
            if ~ truncatedView( il )
                leftTA = totalAttenuations( il );
                break;
            end
        end
        
                ir = iv;
        while ir <  noViews
            ir = ir + 1;
            if ~ truncatedView( ir )
                rightTA = totalAttenuations( ir );
                break;
            end
        end
        
        if leftTA == 0
            totalAttenuations( iv ) = rightTA;
        elseif rightTA == 0
            totalAttenuations( iv ) = leftTA;
        else
            totalAttenuations( iv ) = ( ( ir - iv ) * leftTA + ( iv - il ) * rightTA ) / ( ir - il );
        end
                   
    end
    
end





end
