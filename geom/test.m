% Example of multi threading in a Mex File

% Compile the mex file
mex square.c -v;

% Input array
x = [1 2 3 4 5 6 7 8 9];

% Execute the multi-threaded mex file, which will square all x values
y = square(x);

% Display the squared values
disp(y);

% Demonstrate with larger matrix to show (realtime) percentage information
y = square(rand(100,100));
