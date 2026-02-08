%% Z-factor calculation for cubic EOS


function zfactor = calczfactor(a, b, eos_type)

% Default to PR if eos_type not specified
if nargin < 3 || isempty(eos_type)
    eos_type = 'PR';
end


switch upper(eos_type)
    case 'PR'
        % Z^3 - (1-B)Z^2 + (A - 3B^2 - 2B)Z - (AB - B^2 - B^3) = 0
        c1 = 1;
        c2 = b - 1;
        c3 = a - 3*b^2 - 2*b;
        c4 = -a*b + b^2 + b^3;
        
    case 'SRK'
        % Z^3 - Z^2 + (A - B - B^2)Z - AB = 0
        c1 = 1;
        c2 = -1;
        c3 = a - b - b^2;
        c4 = -a*b;
end


zroots = roots([c1 c2 c3 c4]);

% Choose the real roots (with small tolerance for numerical noise)
tol = 1e-10;
zfactor = [];
for i = 1:3
    if abs(imag(zroots(i))) < tol
        z = real(zroots(i));
        if z > 0  % Z must be positive
            zfactor = cat(1, zfactor, z);
        end
    end
end

zfactor = sort(zfactor);

end