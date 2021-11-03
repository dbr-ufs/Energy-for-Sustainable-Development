syms eta_ref beta T_cel T_re FR_Tau_a FR_UL T_amb Irrad T_fria T_gelada COP_chiller
syms eta_ther eta_el eta_1 eta_2 eta_frio eta_2chiller COP_e

eqn(1) = eta_el == eta_ref*(1-beta*(T_cel-25)); % T_ref=25;
eqn(2) = eta_ther==FR_Tau_a*(1-eta_el)-FR_UL*(T_cel-T_amb)/Irrad;
eqn(3) = eta_1==eta_el+eta_ther;
eqn(4) = eta_2==(1-(T_amb+273)/(T_cel+273))*eta_ther+eta_el;
eqn(5) = eta_frio==eta_ther*COP_chiller;
eqn(6) = eta_2chiller==((T_amb+273)/(T_gelada+273)-1)*eta_frio+eta_el;
%eqn(7) = eta_2==eta_ref; % máxima temperatura para aproveitamento do calor
eqn(7) = eta_2chiller==eta_ref; % busca-se aqui o melhor aproveitamento do recurso solar
%eqn(7) = eta_2chiller==eta_ref-eta_frio/COP_e;  %condição de comparação igual, quando ambos os recursos são necessários 
eqn(8) = beta == 0.005;
eqn(9) = T_fria==30;
eqn(10) = T_gelada==5;
eqn(11) = COP_chiller==0.55;
eqn(12) = T_amb==28;
eqn(13) = Irrad==1000;
eqn(14) = FR_Tau_a==0.574/(1-0.123);
eqn(15) = FR_UL==4.85;
eqn(16) = eta_ref==0.2;
eqn(17) = COP_e == 3;
clear T_max Ir

for k = 1:5
 bet(k) = 0.001+(k-1)*0.001;
 eqn(8) = beta == bet(k);
 for j=1:3
  eta_r(j) = 0.05+(j-1)/1*0.05;        eta_r(3) = 0.2;
  eqn(16) = eta_ref==eta_r(j);
  for i=1:10
    Ir(i) = (i-1)^2/81*10000;
    eqn(13) = Irrad==Ir(i);
    S=vpasolve(eqn);
    if max(S.T_cel)>0 & isreal(max(S.T_cel))
        T_max(i,3*(k-1)+j)=max(S.T_cel);
    else
        T_max(i,3*(k-1)+j) = 0;
    end
  end
 end
end
T_max = double(T_max);
%gera_graf_com_frio(Ir,T_max)
%gera_graf_com_calor(Ir,T_max)
%Tmax_Com_graf(Ir,T_max)
Tmax_Com_graf_frio(Ir,T_max)