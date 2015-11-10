function sinoQuick( sino )

h = figure;
if ndims( sino ) == 3
   
    imdisp( sino(end/2,:,end/4:end/2 ) );
else
    imdisp( sino  );
end

set( h, 'NumberTitle', 'off', 'Name', [inputname(1) ' ' num2str( fix( mod( cputime * 100, 10000) ) )]);
end