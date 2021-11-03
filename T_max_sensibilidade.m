syms eta_ref beta T_cel T_re FR_Tau_a FR_UL T_amb Irrad T_fria T_gelada COP_chiller
syms eta_ther eta_el eta_1 eta_2 eta_frio eta_2chiller COP_e

%%CALOR
% eqn = eta_ref==(1-(T_amb+273)/(T_cel+273))*FR_Tau_a*(1-eta_ref*(1-beta*(T_cel-25)))-FR_UL*(T_cel-T_amb)/Irrad+eta_ref*(1-beta*(T_cel-25));
% S1=solve(eqn,T_cel);
% S2=S1(2);
% S1=S1(1);

%FRIO
eqn = eta_ref==((T_amb+273)/(T_gelada+273)-1)*FR_Tau_a*(1-eta_ref*(1-beta*(T_cel-25)))-FR_UL*(T_cel-T_amb)/Irrad*COP_chiller+eta_ref*(1-beta*(T_cel-25));
S1=solve(eqn,T_cel);

eqns(2) = beta == 0.001;
eqns(3) = T_fria==30;
eqns(4) = T_gelada==5;
eqns(5) = COP_chiller==0.01;%(0.003426*T_cel^3  -0.2315*T_cel^2 + 4.074*T_cel - 11.06) / (T_cel^2 -98.66*T_cel + 2574);%0.55;
eqns(6) = T_amb==28;
eqns(7) = Irrad==1000;
eqns(8) = FR_Tau_a==0.574/(1-0.123);
eqns(9) = FR_UL==4.85;
eqns(10) = eta_ref==0.1;

eqns(1) = T_cel == simplify(eta_ref/S1*diff(S1,eta_ref));
for i = 1:11
    Ir(i) = (i-1)^2/100*10000; %4000+(i-1)/10*6000;%
    eqns(7) = Irrad==Ir(i);        
    S=vpasolve(eqns);
    %if double(S.T_cel) < 60
    Sensib(i,1)=double(S.T_cel);
end

eqns(1) = T_cel == simplify(Irrad/S1*diff(S1,Irrad));
for i = 1:11
    eqns(7) = Irrad==Ir(i);
    S=vpasolve(eqns);
    Sensib(i,2)=double(S.T_cel);
end

eqns(1) = T_cel == simplify(beta/S1*diff(S1,beta));
for i = 1:11
    eqns(7) = Irrad==Ir(i);
    S=vpasolve(eqns);
    Sensib(i,3)=double(S.T_cel);
end

eqns(1) = T_cel == simplify(FR_Tau_a/S1*diff(S1,FR_Tau_a));
for i = 1:11
    eqns(7) = Irrad==Ir(i);
    S=vpasolve(eqns);
    Sensib(i,4)=double(S.T_cel);
end

eqns(1) = T_cel == simplify(FR_UL/S1*diff(S1,FR_UL));
for i = 1:11
    eqns(7) = Irrad==Ir(i);
    S=vpasolve(eqns);
    Sensib(i,5)=double(S.T_cel);
end

eqns(1) = T_cel == simplify(T_amb/S1*diff(S1,T_amb));
for i = 1:11
    eqns(7) = Irrad==Ir(i);
    S=vpasolve(eqns);
    Sensib(i,6)=double(S.T_cel);
end

eqns(1) = T_cel == simplify(T_gelada/S1*diff(S1,T_gelada));
for i = 1:11
    eqns(7) = Irrad==Ir(i);
    S=vpasolve(eqns);
    Sensib(i,7)=double(S.T_cel);
end

eqns(1) = T_cel == simplify(COP_chiller/S1*diff(S1,COP_chiller));
for i = 1:11
    eqns(7) = Irrad==Ir(i);
    S=vpasolve(eqns);
    Sensib(i,8)=double(S.T_cel);
end

%T_max_sensib_graf(Ir,Sensib*100) 
T_max_sensib_com_graf(Ir,Sensib*100)