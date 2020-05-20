function sig = GenSignal(sigOpt)

% switch(sigOpt)
%     case 1
%         sig = @(t) sin(2*t).^5.*(0<=t);
%     case 2
%         sig = @(t) sin(16*t).^5.*(0<=t).*(t<=0.40);
%     case 3
%         sig = @(t) sin(32*t).^5.*(0<=t);
%     case 4
%         sig = @(t) sin(32*t).^5.*(0<=t).*(t<=0.49);
% end

switch(sigOpt)
    case 1
        sigStr = '@(t) sin(2*t).^5.*(0<=t)';
    case 2
        sigStr = '@(t) sin(16*t).^5.*(0<=t).*(t<=0.40)';
    case 3
        sigStr = '@(t) sin(32*t).^5.*(0<=t)';
    case 4
        sigStr = '@(t) sin(32*t).^5.*(0<=t).*(t<=0.49)';
end
    
sig = str2func(sigStr);
fprintf('Generated signal: %s\n', sigStr)

end