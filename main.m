
%% SETUP
addpath('periodic-euler','helper-functions')

% FOM parameters
N       = 200;
dt      = 1e-5;
Tfinal  = 1e-2;
xgrid   = linspace(0,2,N+1)';
xgrid   = xgrid(1:end-1);

% non-dimensionalization scaling constants
ndc.rho = 10;   
ndc.u = 100;
ndc.l = 2;
ndc.t = ndc.l/ndc.u;

%% ACQUIRE TRAINING DATA FROM NONLINEAR FOM
% run FOM to acquire training state data
M = 20;
s_train = cell(M,1);
for i = 1:M
    init = rand_init(xgrid);
    [s,init,time,runtime] = FOM(xgrid,dt,init,Tfinal);
    s_train{i} = [init.s0, s];
end

%% TRANSFORM DATA TO QUADRATIC STATE
% process data to get non-dim transformed state & time-derivative data
w_all = cell(M,1);
r_all = cell(M,1);
for i = 1:M
    % transform state to specific volume representation
    w        = nonlin2quad(s_train{i});     
    
    % non-dimensionalize state so that POD basis is better
    w_nd     = nondim(w,ndc,'spv');
    
    % compute non-dim residual in sp. vol. variables
    r_nd     = (w_nd(:,2:end)-w_nd(:,1:end-1))/(dt/ndc.t);
    
    w_all{i} = w_nd(:,1:end-1);
    r_all{i} = r_nd;
end
W = cat(2,w_all{:});
R = cat(2,r_all{:});

%% FORMULATE AND SOLVE LEARNING PROBLEM
% take SVD to get POD modes
[U,s,~] = svd(W,'econ');

% set POD basis size (same as dimension of model to be learned)
r = 20;
Ur = U(:,1:r);

% project data
What = Ur'*W;
Rhat = Ur'*R;

% get quadratic data w/o Kronecker redundancy
Wsq = get_x_sq(What');

% solve learning problem
F = tikhonov(Rhat',Wsq,0.001)';    % gets non-redundant coefficients of quadratic terms
H = F2H(F);                         % converts coefficients to matricized tensor H

disp(['Transform & Learn model of size r = ',num2str(r),' learned from ',num2str(M),' random training trajectories.'])
%% COMPUTE TRANFORM & LEARN TRAINING ERROR
ndtime = time/ndc.t;
f_rom = @(w) H*kron(w,w);     % define nondim ROM residual

training_err = zeros(M,1);
for i = 1:M
    what0 = Ur'*w_all{i}(:,1);                  % ROM initial condition
    w_rom = forwardEuler(f_rom,what0,ndtime);   % non-dimensional integration
    w_rec = Ur*w_rom;                           % reconstruct high-dim state
    
    % compute reference solution in nondim sp. vol. variables
    w_ref = nondim(nonlin2quad(s_train{i}(:,2:end)),ndc,'spv');
    
    training_err(i) = norm(w_rec - w_ref,'fro')/norm(w_ref,'fro');
end

disp(['Mean relative training error: ',num2str(mean(training_err))])

%% COMPUTE TEST DATA SET
s_test = cell(M,1);
for i = 1:M
    init = rand_init(xgrid);
    [s,init,time,runtime] = FOM(xgrid,dt,init,Tfinal);
    s_test{i} = [init.s0, s];
end

%% COMPUTE TEST ERROR
test_err = zeros(M,1);
for i = 1:M
    % ROM initial condition 
    what0 = Ur'*nondim(nonlin2quad(s_test{i}(:,1)),ndc,'spv');  
     
    w_rom = forwardEuler(f_rom,what0,ndtime);   % non-dimensional integration
    w_rec = Ur*w_rom;                           % reconstruct high-dim state
    
    w_ref = nondim(nonlin2quad(s_test{i}(:,2:end)),ndc,'spv');
    test_err(i) = norm(w_rec - w_ref,'fro')/norm(w_ref,'fro');
end

disp(['Mean relative test error:     ',num2str(mean(test_err))])