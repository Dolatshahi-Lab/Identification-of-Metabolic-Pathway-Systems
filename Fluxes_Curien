function [fluxes,xdot] = Fluxes_Curien(t, y)

% ----------------------------------------------------------------------- %
%
% This function outputs the fluxes and xdots for comparison with the result
% of analyses.
% This model and the parameters are taken from the Curien paper and rewritten
% with MATLAB syntax. In this model adoption, Vcgs was set to zero since it 
% was more than 5 times smaller than the smallest flux and in order to keep
% the degrees of freedom at 2 for better visualization.
%
% ----------------------------------------------------------------------- %

% ------------------------- parameters ---------------------------------- %
AK1 = 0.25;
AK2 = 0.25;
AKHSDHII = 0.63;
AKHSDHI = 0.63;
ASADH = 11.6;
DHDPS1 = 1.6; 
DHDPS2=1.6;
HSK = 4;
TS1 = 7.4;
TD=0.36;
CGS = 0.7;
THA=0;
LKR=0;
AdoMet=20;
Cys=15;
Phosphate=10000;
Val=100;
AK1_kforward_app_exp=5.65;
AK1_kreverse_app_exp=1.57;
AK1_Lys_Ki_app_exp=550;
AK1_AdoMet_Ka_app_exp=3.5;
AK1_h_exp=2;
AK2_kforward_app_exp=3.15;
AK2_kreverse_app_exp=0.88;
AK2_Lys_Ki_app_exp=22;
AK2_h_exp=1.1;
AKI_kforward_app_exp=0.36;
AKI_kreverse_app_exp=0.10;
AKI_Thr_Ki_app_exp=124;
AKI_h_exp=2;
AKII_kforward_app_exp=1.35;
AKII_kreverse_app_exp=0.38;
AKII_Thr_Ki_app_exp=109;
AKII_h_exp=2;
ASADH_kforward_app_exp=0.9;
ASADH_kreverse_app_exp=0.23;
DHDPS1_k_app_exp=1;
DHDPS1_Lys_Ki_app_exp=10;
DHDPS1_h_exp=2;
DHDPS2_k_app_exp=1;
DHDPS2_Lys_Ki_app_exp=33;
DHDPS2_h_exp=2;
HSDHI_kforward_app_exp=0.84;
HSDHI_Thr_relative_residual_activity_app_exp=0.15;
HSDHI_Thr_relative_inhibition_app_exp=0.85;
HSDHI_Thr_Ki_app_exp=400;
HSDHII_kforward_app_exp=0.64;
HSDHII_Thr_relative_residual_activity_app_exp=0.25;
HSDHII_Thr_relative_inhibition_app_exp=0.75;
HSDHII_Thr_Ki_app_exp=8500;
HSK_kcat_app_exp=2.8;
HSK_Hser_app_exp=14;
TS1_kcatmin_exp=0.42;
TS1_AdoMet_kcatmax_exp=3.5;
TS1_AdoMEt_Km_no_AdoMet_exp=250;
TS1_AdoMet_Ka1_exp=73;
TS1_AdoMet_Ka2_exp=0.5;
TS1_AdoMet_Ka3_exp=1.09;
TS1_AdoMet_Ka4_exp=142;
TS1_Phosphate_Ki_exp=1000;
TS1_h_exp=2;
CGS_kcat_exp=30;
CGS_Cys_Km_exp=460;
CGS_Phser_Km_exp=2500;
CGS_Phosphate_Ki_exp=2000;
TD_k_app_exp=0.0124;
TD_Ile_Ki_no_Val_app_exp=30;
TD_Val_Ka1_app_exp=73;
TD_Val_Ka2_app_exp=615;
TD_h_app_exp=3;
Lys_tRNAS_Vmax=0.43;
Lys_tRNAS_Lys_Km=25;
Thr_tRNAS_Vmax=0.43;
Thr_tRNAS_Thr_Km=100;
Ile_tRNAS_Vmax=0.43;
Ile_tRNAS_Ile_Km=20;
THA_kcat_exp=1.7;
THA_Thr_Km_exp=7100;
LKR_kcat_exp=3.1;
LKR_Lys_Km_exp=13000;


% --------------------  Flux formulations --------------------------------% 
Vak1 =AK1*(AK1_kforward_app_exp - AK1_kreverse_app_exp*y(1,:))./ ...
    (1+(y(3,:)/(AK1_Lys_Ki_app_exp/(1+AdoMet/ AK1_AdoMet_Ka_app_exp))).^AK1_h_exp);
Vak2 = AK2*(AK2_kforward_app_exp - AK2_kreverse_app_exp*y(1,:))./ ...
    (1+(y(3,:)/ AK2_Lys_Ki_app_exp).^AK2_h_exp);
VakI =AKHSDHI*(AKI_kforward_app_exp - AKI_kreverse_app_exp*y(1,:)).*1./ ...
    (1+(y(6,:)/ AKI_Thr_Ki_app_exp).^AKI_h_exp);
VakII =AKHSDHII*(AKII_kforward_app_exp - AKII_kreverse_app_exp*y(1,:))./...
    (1+(y(6,:)/ AKII_Thr_Ki_app_exp).^AKII_h_exp);
Vasadh = ASADH*(ASADH_kforward_app_exp*y(1,:)- ASADH_kreverse_app_exp*y(2,:));
Vdhdps1 = DHDPS1* DHDPS1_k_app_exp *y(2,:).*(1./(1+(y(3,:)/DHDPS1_Lys_Ki_app_exp).^DHDPS1_h_exp));
Vdhdps2 = DHDPS2* DHDPS2_k_app_exp *y(2,:).*(1./(1+(y(3,:)/ DHDPS2_Lys_Ki_app_exp).^DHDPS2_h_exp));
VhsdhI = AKHSDHI* HSDHI_kforward_app_exp *y(2,:).*( HSDHI_Thr_relative_residual_activity_app_exp + ...
    HSDHI_Thr_relative_inhibition_app_exp ./(1+y(6,:)/ HSDHI_Thr_Ki_app_exp));
VhsdhII = AKHSDHII* HSDHII_kforward_app_exp *y(2,:).*( HSDHII_Thr_relative_residual_activity_app_exp + ...
    HSDHII_Thr_relative_inhibition_app_exp ./(1+y(6,:)/ HSDHII_Thr_Ki_app_exp));
Vhsk = HSK* HSK_kcat_app_exp *y(4,:)./( HSK_Hser_app_exp +y(4,:));
Vts1=TS1*(TS1_kcatmin_exp + TS1_AdoMet_kcatmax_exp *AdoMet^TS1_h_exp /TS1_AdoMet_Ka1_exp)/ ...
    (1+AdoMet^TS1_h_exp / TS1_AdoMet_Ka1_exp)*y(5,:)./((1+Phosphate/TS1_Phosphate_Ki_exp)* ...
    ( TS1_AdoMEt_Km_no_AdoMet_exp *(1+AdoMet/ TS1_AdoMet_Ka2_exp)/...
    (1+AdoMet/ TS1_AdoMet_Ka3_exp))./(1+AdoMet^TS1_h_exp / TS1_AdoMet_Ka4_exp)+y(5,:));
% Vcgs= CGS*( CGS_kcat_exp /(1+ CGS_Cys_Km_exp /Cys))*y(5,:)./...
%     (( CGS_Phser_Km_exp /(1+ CGS_Cys_Km_exp /Cys)).*(1+Phosphate/ CGS_Phosphate_Ki_exp)+ y(5,:));
Vtd = TD*TD_k_app_exp*y(6,:)./(1+(y(7,:)/( TD_Ile_Ki_no_Val_app_exp + TD_Val_Ka1_app_exp *Val/...
    ( TD_Val_Ka2_app_exp +Val))).^ TD_h_app_exp);
VLys_tRNAS = Lys_tRNAS_Vmax *y(3,:)./(Lys_tRNAS_Lys_Km+y(3,:));
VThr_tRNAS = Thr_tRNAS_Vmax *y(6,:)./(Thr_tRNAS_Thr_Km+y(6,:));
VIle_tRNAS =Ile_tRNAS_Vmax *y(7,:)./(Ile_tRNAS_Ile_Km+y(7,:));
Vlkr= LKR*LKR_kcat_exp*y(3,:)./(LKR_Lys_Km_exp+y(3,:));

% ------------------ fluxes are the main outputs of this function ------- % 
fluxes = [Vak1 + Vak2 + VakI + VakII;Vasadh;VhsdhI+VhsdhII;Vdhdps1+Vdhdps2;...
    Vlkr+VLys_tRNAS;Vhsk;Vts1;VThr_tRNAS;Vtd;VIle_tRNAS];

% --- xdots are the secondary outputs of this function in case needed for 
% trouble shooting ------------------------------------------------------ % 
xdot = zeros(size(y));

% dydt(1) = vAK - vASADH ;
xdot(1,:) = Vak1 + Vak2 + VakI + VakII - Vasadh;
% dydt(2) = vASADH - vDHDPS - vHSDH ;
xdot(2,:) = Vasadh - Vdhdps1 -Vdhdps2- VhsdhI -VhsdhII ;
% dydt(3) = vDHDPS - vLKR ;
xdot(3,:) = Vdhdps1 + Vdhdps2 - Vlkr - VLys_tRNAS;
% dydt(4) = vHSDH - vHSK ;
xdot(4,:) = VhsdhI + VhsdhII -Vhsk;
% dydt(5) = vHSK - vTS1 - Vcgs;
xdot(5,:) = Vhsk  - Vts1;
% dydt(6) = vTS1 - vTD - vThr_tRNAsth;
xdot(6,:) = Vts1 - Vtd - VThr_tRNAS;
% dydt(7) = Vtd-VIle_tRNAS;
xdot(7,:) = Vtd-VIle_tRNAS;
% dydt(8) = vThr_tRNAsth;
xdot(8,:) = VThr_tRNAS;
end
