%% General settings

param.tau = 1.8;
param.T1b = 1.65;
param.T1t = 1.445;
param.noise = 0.001545;
param.M0B = 1;
param.lamda = 0.9;
param.f = 50/6000;
param.alpha = 0.85;
stepSize = 0.001;
BATDist = 0.2:stepSize:2.3;
distWeight = [linspace(0.5,1,length(BATDist(1):stepSize:0.5-stepSize)), ones(1,length(0.5:stepSize:2)), linspace(1,0.5,length(2+stepSize:stepSize:BATDist(end)))];
distWeight = distWeight((length(distWeight)-length(BATDist)+1):end);
scanTime = 300;
allSlice = 1;
slicedt = 0;
A = [1, 0 ;0, 0];
lims.PLDStep = 0.025;
lims.PLDLB = 0.075;
lims.tauStep = 0.025;
lims.tauLB = 0.1;
lims.tauUB = 1.8;
tReadout = 0.5;

%% Seq_singleLD and Seq_multiLD

rng('default')
param.filename = 'var_multi_pCASL';
lims.PLDUB = BATDist(end)+lims.PLDStep;
nPLD = 10;
[bestPLD,bestTau,bestminVariance] = OED_PCASL_SeqsingleLD_LOptimal(param,BATDist,distWeight,scanTime,A,allSlice,nPLD,lims,slicedt);

% [bestPLD,bestTau,bestminVariance] = OED_PCASL_SeqmultiLD_LOptimal(param,BATDist,distWeight,scanTime,A,allSlice,nPLD,lims,slicedt);

%% Time-encoded (fixed or T1-adjusted duration or free-lunch)

param.filename = 'var_te_pCASL';
lims.PLDUB = 1;
nPLD = 7;
param.num_enc = nPLD;
param.multiPLD = 1;

[bestPLD,bestTau,bestminVariance] = OED_PCASL_Hadfixed_LOptimal(param,BATDist,distWeight,scanTime,A,allSlice,nPLD,lims,slicedt);

% [bestPLD,bestTau,bestminVariance] = OED_PCASL_HadT1adj_LOptimal(param,BATDist,distWeight,scanTime,A,allSlice,nPLD,lims,slicedt);

% [bestPLD,bestTau,bestminVariance] = OED_PCASL_Hadfreelunch_fixed_LOptimal(param,BATDist,distWeight,scanTime,A,allSlice,nPLD,lims,slicedt);

% [bestPLD,bestTau,bestminVariance] = OED_PCASL_Hadfreelunch_T1adj_LOptimal(param,BATDist,distWeight,scanTime,A,allSlice,nPLD,lims,slicedt);

% rng('default')
% [bestPLD,bestTau,bestminVariance] = OED_PCASL_Hadvariable_LOptimal(param,BATDist,distWeight,scanTime,A,allSlice,nPLD,lims,slicedt);


%% Hybrid (fixed or T1-adjusted duration)

param.filename = 'var_te_pCASL_nPLD';
lims.PLDUB = 1;
nPLD = 3;
param.num_enc = nPLD;
param.multiPLD = 3;

% [bestPLD,bestTau,bestminVariance] = OED_PCASL_Hybridfixed_LOptimal(param,BATDist,distWeight,scanTime,A,allSlice,nPLD,lims,slicedt);

[bestPLD,bestTau,bestminVariance] = OED_PCASL_HybridT1adj_LOptimal(param,BATDist,distWeight,scanTime,A,allSlice,nPLD,lims,slicedt);

% rng('default')
% [bestPLD,bestTau,bestminVariance] = OED_PCASL_Hybridvariable_LOptimal(param,BATDist,distWeight,scanTime,A,allSlice,nPLD,lims,slicedt);

