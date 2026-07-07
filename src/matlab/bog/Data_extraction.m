% Data_extraction.m
% Post-processes output from LH2Simulate.
%
% Configuration:
%   (ST) or 1 = Storage tank at liquefaction terminal  (FEEDING vessel)
%   (ET) or 2 = LH2 trailer                           (RECEIVING vessel)


%nominal.rho_L1 = -0.1328061400 * (nominal.uL1(:,LH2Model.nL1))./1000 + 70.8804573705; % correlation from REFPROP V9.1
nominal.rho_L1 = -5.12074746E-07*((nominal.uL1(:,LH2Model.nL1))./1000).^3 - 1.56628367E-05*((nominal.uL1(:,LH2Model.nL1))./1000).^2 - 1.18436797E-01*((nominal.uL1(:,LH2Model.nL1))./1000) + 7.06218354E+01;
nominal.VL1 = nominal.mL1./nominal.rho_L1;
nominal.Vullage1 = LH2Model.VTotal1-nominal.VL1;
nominal.rhov1 = nominal.mv1./nominal.Vullage1;

%nominal.rho_L2 = -0.1328061400 * (nominal.uL2(:,LH2Model.nL2))./1000 + 70.8804573705; % correlation from REFPROP V9.1
nominal.rho_L2 = -5.12074746E-07*((nominal.uL2(:,LH2Model.nL2))./1000).^3 - 1.56628367E-05*((nominal.uL2(:,LH2Model.nL2))./1000).^2 - 1.18436797E-01*((nominal.uL2(:,LH2Model.nL2))./1000) + 7.06218354E+01; % correlation from REFPROP v9.1
  
nominal.VL2 = nominal.mL2./nominal.rho_L2;

% Compute liquid height for horizontal ET using cylVToH
nominal.hL2 = zeros(size(nominal.VL2));
for zz = 1:length(nominal.VL2)
    nominal.hL2(zz) = cylVToH(nominal.VL2(zz), LH2Model.R2, LH2Model.Lcyl2);
end

nominal.Vullage2= LH2Model.VTotal2-nominal.VL2;
nominal.rhov2 = nominal.mv2./nominal.Vullage2;
nominal.pcthL2 = nominal.VL2 ./ LH2Model.VTotal2; % volume-based fill fraction
    
for z=1:length(nominal.rhov1);
    
        nominal.pv1(z,:)=vaporpressure_bis(nominal.uv1(z,end), nominal.rhov1(z));
        for ii = 1:LH2Model.nV1
            if nominal.uv1(z,ii)<0*(-108.8/(2.0159/1000))
                nominal.uv1(z,ii)= 0;
            end

            nominal.Tv1(z,ii) = refpropm('T','P',nominal.pv1(z)/1000,'U',nominal.uv1(z,ii),'PARAHYD');
        end
        for ii = 1:LH2Model.nL1
            nominal.TL1(z,ii) = -0.0002041552*(nominal.uL1(z,ii)/1000)^2 + 0.1010598604*nominal.uL1(z,ii)/1000 + 20.3899281428;
        end
        nominal.Jvalve111(z,:)=gasFlow(LH2Model.S_valve1,LH2Model.gamma_,nominal.rhov1(z),nominal.pv1(z),LH2Model.p_atm);
        nominal.hL1(z) = cylVToH(nominal.VL1(z),LH2Model.R1,LH2Model.Lcyl);
  
       
        nominal.pv2(z,:)=vaporpressure_bis(nominal.uv2(z,end),nominal.rhov2(z));
        nominal.Jvalve222(z,:)=nominal.ETTTVenstate(z)*gasFlow(LH2Model.S_valve2,LH2Model.gamma_,nominal.rhov2(z),nominal.pv2(z),LH2Model.p_atm);
        for ii = 1:LH2Model.nV2
            if nominal.uv2(z,ii)<0*(-108.8/(2.0159/1000))
                nominal.uv2(z,ii)= 0;
            end
            nominal.Tv2(z,ii) = refpropm('T','P',nominal.pv2(z)/1000,'U',nominal.uv2(z,ii),'PARAHYD');
        end
        for ii = 1:LH2Model.nL2
            nominal.TL2(z,ii) = -0.0002041552*(nominal.uL2(z,ii)/1000)^2 + 0.1010598604*nominal.uL2(z,ii)/1000 + 20.3899281428;
        end
end
    
    nominal.pL1 = nominal.rho_L1.*LH2Model.g.*(nominal.hL1)';
    nominal.pTotal1 = nominal.pv1+nominal.pL1;
    
    % CORRECTED: Interface area for horizontal ST (chord * length),
    % consistent with the computation inside LH2Simulate.m.
    nominal.S1 = zeros(size(nominal.hL1));
    for zz = 1:length(nominal.hL1)
        hh = nominal.hL1(zz);
        if hh > LH2Model.R1
            dd = hh - LH2Model.R1;
        else
            dd = LH2Model.R1 - hh;
        end
        dd = min(dd, LH2Model.R1 * 0.999); % clamp for edge cases
        cc = 2 * LH2Model.R1 * sqrt(1 - (dd/LH2Model.R1)^2);
        nominal.S1(zz) = cc * LH2Model.Lcyl;
    end
      
    nominal.pL2 = nominal.rho_L2.*LH2Model.g.*(nominal.hL2);
    nominal.pTotal2 = nominal.pv2+nominal.pL2;
    

    nominal.Boiloff_ET=zeros(length(nominal.t),1);
    for ii = 2:length(nominal.t)
        nominal.Boiloff_ET(ii) = nominal.Boiloff_ET(ii-1) + (nominal.t(ii)-nominal.t(ii-1))* nominal.Jvalve222(ii,:);
    end
    
    nominal.Jv10 = nominal.Jboil -nominal.AAA - nominal.Jcd1;
    nominal.Jv20 = -nominal.Jvalve222 - nominal.Jcd2;
    display('Data extraction done');
