%% CALCULATE FUGACITY COEFFICIENT FOR MULTI-COMPONENT SYSTEMS
% -------------------------------------------------------------------------
% The Definition of Variables.
% zfactor  : compressibility factor
% A        : component dimensionless attraction parameters
% B        : component dimensionless covolume parameters  
% Amix     : mixture dimensionless attraction parameter
% Bmix     : mixture dimensionless covolume parameter
% Amix2    : Amix2(i) = sum_j(x_j * A_ij)
% eos_type : 'PR' (default) or 'SRK'
% -------------------------------------------------------------------------

function fugcoef = calcfugcoef_multicomp(zfactor, A, B, Amix, Bmix, Amix2, eos_type)

% Default to PR if eos_type not specified
if nargin < 7 || isempty(eos_type)
    eos_type = 'PR';
end

ncomp = size(A,1);
fugcoef = zeros(ncomp, 1);

% Ensure zfactor is scalar
if length(zfactor) > 1
    zfactor = zfactor(1);  % Take first root if multiple
end

for i = 1:ncomp
    
    if zfactor < Bmix
        error('Z-factor (%.4f) must be larger than Bmix (%.4f).\n', zfactor, Bmix);
    end
    
    switch upper(eos_type)
        case 'PR'
            % Peng-Robinson fugacity coefficient
            c0 = 2*sqrt(2);
            c1 = 1 + sqrt(2);
            c2 = 1 - sqrt(2);
            
            lnphi = B(i)/Bmix*(zfactor - 1) - log(zfactor - Bmix) ...
                - Amix/(c0*Bmix) * (2*Amix2(i)/Amix - B(i)/Bmix) ...
                * log((zfactor + c1*Bmix)/(zfactor + c2*Bmix));
            
        case 'SRK'
            % Soave-Redlich-Kwong fugacity coefficient
            lnphi = B(i)/Bmix*(zfactor - 1) - log(zfactor - Bmix) ...
                - Amix/Bmix * (2*Amix2(i)/Amix - B(i)/Bmix) ...
                * log(1 + Bmix/zfactor);
            
        otherwise
            error('Unknown EOS type: %s. Use ''PR'' or ''SRK''.', eos_type);
    end
    
    fugcoef(i) = exp(lnphi);
    
end

end