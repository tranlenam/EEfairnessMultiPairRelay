clc
rng(1)
%% Simulation Parameters
noS = 3; % No of source nodes
noD = 3; % No of destination nodes
noR = 8; % No of relay nodes
No = -90 -30; % Noise density (dBW)
BW = 250e3; % Bandwidth
T_s = 10e-3; % 10ms
P_max = 10^(33/10)/1000 ; % Users transmit power 
Q0 = 01; % required User QoS
P_cir1 =(0.1+1.0)/1000; % Users' circuit power P'
P_cir2 =(0.1+1.0)/1000; % Users' circuit power P''
P_Rcir =1/1000; % Relay circuit power

% EH circuit specifications
P_dc_EH = 24/1000; 
c = 150;
d = 0.014;
myalpha_0 = exp(c*d);
mybeta_0 = 1/(1+exp(c*d));
mybeta_1 = P_dc_EH/(1-mybeta_0);

% PA parameter
max_PA_eff = 0.35;
power_PA_max = 10^(33/10)/1000;
eps_R = sqrt( power_PA_max)/max_PA_eff;
eps_S = sqrt(P_max)/max_PA_eff;
% SP power coefficient
p_sp = 50e-12*BW; 

% channel generation
load ../data/chantest.mat;
%% scale
scaleG = 100;
scaleF = 100;
G=G*scaleG;
F=F*scaleF;
sigma_R = scaleF*sqrt(10^((No)/10));
sigma_D = scaleF*scaleG*sqrt(10^((No)/10));

alpha_EH = P_dc_EH/(1-mybeta_0);

%% calculating effective channel
H_kj = cell(noS,noS);
for k=1:noS
   h_kj = transpose(G(:,k).*F(:,1:noS));
   for j=1:noS
        H_kj{k,j} = h_kj(j,:)'*h_kj(j,:);
   end
end
R_temp = G'.*transpose(G);
R = cell(noD,1);
for k=1:noD
    R{k} = sigma_R^2*diag(R_temp(k,:));
end
Q = F.*transpose(F');

%% Initialization
% generate a random initial point
mytau_i = 0.3;
txpower_S_i = P_max/2*ones(1,noS);
beamformer_i = sqrt(P_dc_EH/2)*(randn(noR,1)+1i*randn(noR,1));
r_i = zeros(noD,1);
z_i = zeros(noD,1);
s_i = zeros(noD,1);
v_i = zeros(noD,1);
phi_i = zeros(noD,1);
for k=1:noD
    interf_power = 0;
    for j=1:noD
        if j~=k
            interf_power = interf_power + txpower_S_i(j)*beamformer_i'*H_kj{k,j}*beamformer_i ;
        end
    end
    
    noise_power = sigma_D^2 + beamformer_i'*R{k}*beamformer_i;
    SINR = real(txpower_S_i(k)*beamformer_i'*H_kj{k,k}*beamformer_i )/real(interf_power+noise_power);
    r_i(k) = log(1+SINR);
    z_i(k) = sqrt(mytau_i/r_i(k));
    s_i(k) = real(interf_power+noise_power);
    v_i(k) = SINR;
    phi_i(k) = s_i(k)/v_i(k);

end

mytau_i =  (1+mytau_i)/(1-mytau_i);
p_i = 1./txpower_S_i';
u_i = zeros(noR,1);
t_i = zeros(noR,1);
constraints_violation1_i = zeros(noR,1);
for l=1:noR
    u_i(l) = sqrt(real((beamformer_i(l)'*(txpower_S_i*(abs(F(l,:)).^2)')*beamformer_i(l))));
    power_EH =  txpower_S_i*(abs(F(l,:)/scaleF).^2)';
    t_i(l) = (mytau_i-1)/(1+myalpha_0*exp(-c*power_EH));
    constraints_violation1_i(l) = ((eps_R)*norm([u_i(l),beamformer_i(l)*sigma_R])- scaleF*( t_i(l)*mybeta_1 - (mybeta_0*mybeta_1+P_Rcir)*mytau_i + mybeta_0*mybeta_1-P_Rcir ));
end
constraints_violation2_i = 10*r_i;

myeta_i = max(p_sp+P_cir2./r_i+eps_S*z_i.^2./sqrt(p_i)+P_cir1*z_i.^2);

slack1 = 1;
slack2 = 10000;
ops = sdpsettings('solver','mosek','verbose',0,'debug',1);
% declare optimization variables
mytau = sdpvar(1);
myeta = sdpvar(1);
p = sdpvar(noS,1);
beamformer = sdpvar(noR,1,'full','complex');
z = sdpvar(noD,1);
s = sdpvar(noD,1);
r = sdpvar(noD,1);
v = sdpvar(noD,1);
t = sdpvar(noR,1);
u = sdpvar(noR,1);
% additional variables (for the ease of simulation with Yalmip)
var1 = sdpvar(noR,1);
var2 = sdpvar(noS,1);
var3 = sdpvar(noS,1);
var4 = sdpvar(noS,noD);
var5 = sdpvar(noD,1);
var6 = sdpvar(noR,noD);
var7 = sdpvar(noR,1);
var8 = sdpvar(noD,1);
var9 = sdpvar(noD,1);
var10 = sdpvar(noR,1);
constraints_violation1 = sdpvar(noR,1);
constraints_violation2 = sdpvar(noD,1);

% fixed constraints
C = [];
C = [C; myeta >= 0];
C = [C; mytau >= 1];
C = [C; p >= 1/P_max]; % (43c)
C = [C; z >= 0];
C = [C; t >= 0];
C = [C; u >= 0];
C = [C; r >= 0];
C = [C; s >= 0];
MAXITER = 50;
obj_augmented_seq = zeros(MAXITER,1);
obj_true_seq = zeros(MAXITER,1);
for iIter = 1:MAXITER
    C1=[]; % the constraints that change with iteration
    if nnz( (constraints_violation1_i <= 1e-6)) == noR
        C1 = [C1; constraints_violation1 <= 0]; 
        slack1 = 0; % slack2~=0 implies the finding feasible point process
    end
    
    if nnz( (constraints_violation2_i <= 1e-6) )== noD
        C1 = [C1; constraints_violation2 <= 0];
        slack2 = 0; % slack1~=0 : the finding feasible point process; 
    end
   
    obj_augmented =  myeta +  slack1*sum(constraints_violation2) + slack2*sum(constraints_violation1);
    
    % dyconstraint
        
    C1 = [C1; constraints_violation1 >= 0];
    C1 = [C1; constraints_violation2 >= 0];   
    for k=1:noD 
        % (42)
        C1 = [C1; cone([sqrt(P_cir1)*z(k) (myeta-var8(k)-eps_S*var1(k)-p_sp-1)/2],(myeta-var8(k)-eps_S*var1(k)-p_sp+1)/2)];
        C1 = [C1; cone([z(k) (var1(k)-var2(k))/2],(var1(k)+var2(k))/2)];
        C1 = [C1; cone([var2(k) (p(k)-1)/2],(p(k)+1)/2)];
        C1 = [C1; cone([sqrt(P_cir2) (r(k)-var8(k))/2],(r(k)+var8(k))/2)];
            
        % (44e)        
        C1 = [C1; var3(k) <= 2*(z_i(k)/mytau_i)*z(k) - mytau*(z_i(k)/mytau_i)^2];
        C1 = [C1; cone([1 (var3(k)-r(k))/2],(var3(k)+r(k))/2)];
        
        % (43b)
        C1 = [C1; r(k) + constraints_violation2(k) >= (1+mytau)*Q0];
        
        % (44b)
        C1 = [C1; var5(k) == 2*real(beamformer_i'*H_kj{k,k}'*beamformer)/p_i(k) -  p(k)*real(beamformer_i'*H_kj{k,k}'*beamformer_i)/(p_i(k)^2) ];                       
        C1 = [C1; cone([s(k)/sqrt(2*phi_i(k)), sqrt(phi_i(k)/2)*v(k), (var5(k)-1)/2 ],(var5(k)+1)/2)];
        
        % (35)
        for j=1:noD
            if j~=k
                C1 = [C1; cone([(var4(k,j)-p(j))/2; (H_kj{k,j})^(1/2)*beamformer ],(var4(k,j)+p(j))/2)];
            else 
                C1 = [C1; var4(k,j)==0];
            end
        end
        C1 = [C1; cone([(R{k}')^(0.5)*beamformer; sigma_D; (s(k) -sum(var4(k,:))-1)/2],(s(k) -sum(var4(k,:))+1)/2)];

        % (34)
        C1 = [C1; cone([var9(k), v(k)/2],(v(k)+2)/2)];
        C1 = [C1; cone([(1/sqrt(1+v_i(k))*(log(1+v_i(k))+2-r(k))-var9(k)/2)/2, 1], (1/sqrt(1+v_i(k))*(log(1+v_i(k))+2-r(k))+var9(k)/2)/2)];  
    end
    
    for l=1:noR        
        % (44c) 
        for k=1:noS
              C1 = [C1; cone([abs(F(l,k))*beamformer(l),(p(k)-var6(l,k))/2],(p(k)+var6(l,k))/2)];
        end
        C1 = [C1; sum(var6(l,:))<= 2*u(l)*u_i(l)-u_i(l)^2]; 
        
        % (43d)
        C1 = [C1; cone([u(l),beamformer(l)*sigma_R],sqrt(power_PA_max)*scaleF)];
        
        % (44d)
        C1 = [C1; cone([var10(l), (mytau-t(l)-2)/2],(mytau-t(l))/2)];
        C1 = [C1; cone([(1/sqrt(mytau_i-t_i(l)-1)*(log(mytau_i-t_i(l)-1)+2-var7(l))-var10(l)/2)/2, 1],...
                                                      (1/sqrt(mytau_i-t_i(l)-1)*(log(mytau_i-t_i(l)-1)+2-var7(l))+var10(l)/2)/2)];
                                                  
        C1 = [C1; var7(l) == log(myalpha_0*t_i(l))+ t(l)/t_i(l) - 1 - c*sum((2./p_i-p./(p_i.^2)).*(abs(F(l,:)'/scaleF).^2))];   
            
        % (39)
        C1 = [C1; cone([u(l),beamformer(l)*sigma_R],scaleF/(eps_R)*( t(l)*mybeta_1 - (mybeta_0*mybeta_1+P_Rcir)*mytau + mybeta_0*mybeta_1-P_Rcir + constraints_violation1(l)))];
        
    end
        

    % solution
    solution = optimize([C,C1],obj_augmented,ops);
    
    if (solution.problem == 0) ||  (solution.problem == 3) % if successfully solved or numberical issues occur, continue
        %% Update variables to form problem in the next iteration
        obj_augmented_seq(iIter) = double(obj_augmented);
        obj_true_seq(iIter) = double(myeta);
        u_i = double(u);
        t_i = double(t);
        myeta_i = double(myeta);
        p_i = double(p);
        s_i = double(s);
        z_i = double(z);
        r_i = double(r);
        mytau_i = double(mytau);
        v_i = double(v);
        phi_i = s_i./v_i;
        beamformer_i = double(beamformer);
        constraints_violation1_i = double(constraints_violation1);
        constraints_violation2_i = double(constraints_violation2);
    else % something happens when solving the convex problem
        break
    end

end
plot(1./obj_augmented_seq)
ylabel("Energy Efficiency (nat/s/J)")
xlabel("Iteration index")
saveas(gcf, '../../results/ConvergencePlot.png')