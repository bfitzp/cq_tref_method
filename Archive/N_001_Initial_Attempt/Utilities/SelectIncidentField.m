 function [fn, z_0] = SelectIncidentField(k,opt)
switch(opt{1})
    case 'ps'
        switch(opt{2})
            case 0
                z_0 = 0.0 + 0.0i;
            case 1
                z_0 = 0.5 - 0.5i;
            case 2
                z_0 = -0.5 - 0.5i;
            case 3
                z_0 = 0.95 + 0.95i;
        end
        fn = @(z) 1i/4*besselh(0,k*abs(z-z_0));
    case 'fn'
        switch(opt{2})
            case 1
                fn = @(z) real(z).^2;
        end
end

end