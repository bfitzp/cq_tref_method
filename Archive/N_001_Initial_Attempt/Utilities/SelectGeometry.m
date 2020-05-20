function [w, w_bndBox, cornerInds] = SelectGeometry(opts)
switch(opts)
    case 1
        w = [...
            -1-1i,...
            -1+1i,...
            1+1i,...
            1-1i,...
            ];
        cornerInds = [1;2;3;4];
    case 2
        w = [...
            -1-1i,...
            -0.8+0.3i,...
            -0.1+0.8i,...
            1+1i,...
            1-1i,...
            ];
        cornerInds = [1;2;3;4;5];
    case 3
        w = [...
            -1-1i,...
            -1+0i,...
            0+0i,...
            0+1i,...
            1+1i,...
            1-1i,...
            ];
        cornerInds = [1;2;3;4;5;6];
    case 4
        w = [...
            -1-1i,...
            -1+0i,...
            0.0+0.0i,...
            0+1i,...
            1+1i,...
            1+0i,...
            1-1i,...
            0-1i,...
            ];
        cornerInds = [1;2;3;4;5;7];
    case 5
        w = [ ...
            2 ...
            2+1i ...
            1+1i ...
            1+2i ...
            2i ...
            1i ...
            0 ...
            1];
        cornerInds = [1;2;3;4;5;7];
end

w_bndBox = [min(real(w)), max(real(w)), min(imag(w)), max(imag(w))];    

smallShift = 0; % 1e-12;
w_bndBox(1) = w_bndBox(1) + smallShift; 
w_bndBox(3) = w_bndBox(3) + smallShift;


        
end