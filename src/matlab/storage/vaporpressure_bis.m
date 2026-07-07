function pv = vaporpressure_bis(uv, rhov)
% Standalone version of the vapor pressure function that takes
% internal energy and density as inputs (used by Data_extraction.m).
% Mirrors the nested vaporpressure function in LH2Simulate.m.
%
% Inputs:
%   uv   - internal energy of vapor [J/kg]
%   rhov - density of vapor [kg/m^3]
% Output:
%   pv   - vapor pressure [Pa]

    rhov = max(rhov, 0.01);
    R_H2 = 4124;  % [J/kg/K] specific gas constant for H2
    cv_H2 = 6490; % [J/kg/K] approximate Cv for para-H2 vapor

    quality = 2; % default: assume superheated vapor
    try
        quality = refpropm('q','D',rhov,'U',uv,'PARAHYD');
    catch
        uv_t = fix(100*uv)/100;
        try
            quality = refpropm('q','D',rhov,'U',uv_t,'PARAHYD');
        catch
            rhov_t = max(fix(100*rhov)/100, 0.01);
            try
                quality = refpropm('q','D',rhov_t,'U',uv_t,'PARAHYD');
                rhov = rhov_t;
            catch
                quality = 2;
            end
        end
    end
    
    if quality < 1 && quality > 0
        % Two-phase region
        try
            temp = refpropm('T','D',rhov,'U',uv,'PARAHYD');
            pv = refpropm('P','T',temp,'Q',1,'PARAHYD') * 1e3;
        catch
            T_est = max(uv/cv_H2, 14);
            pv = rhov * R_H2 * T_est;
        end
    else
        % Single-phase (supercritical or superheated)
        try
            pv = refpropm('P','D',rhov,'U',uv,'PARAHYD') * 1e3;
        catch
            T_est = max(uv/cv_H2, 14);
            pv = rhov * R_H2 * T_est;
        end
    end
end
