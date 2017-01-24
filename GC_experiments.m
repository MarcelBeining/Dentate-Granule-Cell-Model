%% comments

%% initialize trees and mechanisms
params = [];

%********* folders 
targetfolder_results = 'C:\Users\beiningm\Dropbox\PhD Deller\Unterlagen\Papers\Paper 3 - P5 compartmental model\Figures\raw';%'D:/GCModel/results';  % the folder where the graphs and pictures are saved
targetfolder_data = 'C:\GCModel\simdata';  % the folder where simulated data is saved (that you do not need to simulate it again if you just want to plot something)
% these params are t2n specific
params.nrnmech = 'nrnmech_win1.dll';   % the name of the compiled mech file (use rename_nrnmech.m to rename your nrnmech.dll file if you use the same lib_mech folder on different machines and thus need different dll files)
params.neuronpath = 'C:\nrn7.4\bin\nrniv.exe';   % the path to your NEURON exe
params.path = 'C:\Users\beiningm\Dropbox\PhD Deller\Modeling\P5 GC Model\aProject Mongiat';  % your main folder of the model
params.morphfolder = 'morphos/NEURON2_hocs';   % folder relative to "path" containing the hoc morphology files
params.exchfolder = 't2nexchange';  % folder name which will be created and is used to exchange data between Matlab and NEURON
% ***********

params.celsius = 24;   % temperature
params.prerun = 400;   % large-dt prerun to let the system equilibrate
params.v_init = -90;  % initial membrane voltage
params.dt = 1;       % standard time step (is automatically changed to a smaller dt in most simulations below)
params.nseg = 'dlambda';  % number of segments, can be constant or 'd_lambda' to adjust it according to the d-lambda rule
params.openNeuron = 0;   % make it 1 to open each NEURON instance (is suppressed if t2n is run with the -q argument)

ostruct = struct('plot','auto','show',3,'legend',0,'marker','o','sem',1,'FontSize',10,'FontType','Arial','barwidth',0.3,'figurewidth',8,'figureheight',5,'LineWidth',1,'grid',0,'priority','plot','ticklength',0.015);  % some options that are used when something is plotted
ostruct.usecol = 1;  % 0 = pseudorandom colors for each simulated cell in graph, 1 = green or blue grading



ostruct.vmodel = 1; % 0 = passive model, > 0 = active model, everything else (e.g. NaN) = old AH99 model
ostruct.changeAHion = 0;  % only important when using the AH99 model. Boolean to decide if standard AH99 ion reversal potentials are used (0) or if they are adjusted to the experiments (1)



ostruct.scalespines = 1;  % scaling of g_pas and cm to implicitly model spines
ostruct.adjustloads = 0;  % the Hay et al 2013 implementation of adjust dendritic loads to reduce variability between cells (not used in publication)
ostruct.noise = 0;       % add noise to the membrane voltage by injecting gaussian noise current to the soma (not working with variable dt / cvode)
%
ostruct.reducecells = 0;  % reduce number of cells for faster simulation (e.g. for testing)
ostruct.usemorph = 1;  % 1 = all SH07, 2= synth mouseMat, 3= synth mouseYoung 4= Beining AAV, 5 = synth ratOld 6= synth ratYoung 7 = Claiborne,
ostruct.newborn = 0;  % 0 = adult GC model, 1 = young abGC model


if ostruct.usemorph >= 4
    ostruct.ratadjust = 1;  % adjust Kir channel in rats
    params.LJP = 0; % no LJP with rat
else
    ostruct.ratadjust = 0;
    params.LJP = 12.1; % liquid junction potential with Mongiat solutions..needs to be substracted from voltage commands and measured voltages
end
ostruct.bablock = 0;  % used in some experiments to distinguish between control (0) and barium application (1) and the corresponding data that should be loaded or simulation that should be done

%*****************************
if ostruct.newborn
%     % this is the good FI,deep fAHP one...
%     ostruct.channelblock = {'Kir21','Kv42','na8st','BK','Cav13','Kv21'};%,'Kv21','Kv42','Kv14','Kv34','na8st','Cav13'};%'Cav22'};     %{'Kv42','Kir21','Kv14','Kv34','pas','Kv21','na8st','Kv723'};%,'na8st','Kv723','na8st','Kv21'};%'except','Kir21'};%{'Kir','Kv42','HCN'};%'Kir','Kv42','Kv14','HCN'}; % Kir , SK etc
%     ostruct.blockamount = [73,50,25,60,50,75];%68 kir %,100,80,80,50,55,100];%  Kv42 50 Kv21 70
%      ostruct.specify = {'','','','','',''};

     %      % this is the used version
    ostruct.channelblock = {'Kir21','Kv42','na8st','BK','Cav13','Kv21','Kv723','BK'};%,'Cav22','Cav12','Cav32'};%,'Kv21','Kv42','Kv14','Kv34','na8st','Cav13'};%'Cav22'};     %{'Kv42','Kir21','Kv14','Kv34','pas','Kv21','na8st','Kv723'};%,'na8st','Kv723','na8st','Kv21'};%'except','Kir21'};%{'Kir','Kv42','HCN'};%'Kir','Kv42','Kv14','HCN'}; % Kir , SK etc
    ostruct.blockamount = [73,50,25,40,50,50,50,100];%,100,100,100];%68 kir %,100,80,80,50,55,100];%  Kv42 50 Kv21 70
     ostruct.specify = {'','','','gakbar','','','','gabkbar'};
    
    ostruct.scalespines = 0.3;  % means g_pas and cm are only scaled by 30% of the original spine densities due to reduced spine density in young abGCs
else
    ostruct.channelblock = {};
    ostruct.blockamount = [];
end

if ~exist(targetfolder_data,'file')
    mkdir(targetfolder_data)
end
if ~exist(targetfolder_results,'file')
    mkdir(targetfolder_results)
end
[tree,params,neuron,treename] = GC_initModel(params,ostruct);  % initialize the model by loading the morphologies and setting the biophysical parameters

neuron_orig = neuron;
params_orig = params;

%% plot trees, Figure 1
% ostruct.show = 2;
% ostruct.savename = 'Fig1_Trees';
% plotmytrees(tree,targetfolder_results,[],ostruct) % '-s'

%% Mongiat IV + Ba, simulate the I-V curve with and without blocking Kir channels with Barium, Figure 2 & 6
neuron = neuron_orig;
ostruct.bablock = 0;
ostruct.holding_voltage = -80; % mV
ostruct.subtract_hv = 1; % boolean subtract holding voltage current
ostruct.show = 1:2; % 1 = only exp data, 2 = only model data, 3 = both
ostruct.variabledt = 1;  % 1 = variable dt, 0 = constant dt
ostruct.coarse = 0;  % 1 = only simulate some voltage steps to make it faster
ostruct.extract_kir = 0; 
ostruct.single = 0;    % 1 = show I-V of each cell
ostruct.figureheight = 4;
ostruct.figurewidth = 6;
%
aGC_IV(neuron,tree,params,targetfolder_data,ostruct);
%
if ostruct.newborn
    ostruct.dataset = 2.28;
    if ~all(ostruct.show==1)
        ostruct.savename = sprintf('Fig6-IV+Ba_young-%s',neuron.experiment);
    else
        ostruct.savename = sprintf('Fig6-IV+Ba_young-%s_onlyexp',neuron.experiment);
    end
else
    ostruct.dataset =3;  % 1 = old mature dataset, 2 = old young dataset, 3 = new mature dataset, 4 = new young dataset, 5 = new mature BaCl dataset, 6 = new young BaCl dataset
    if ~all(ostruct.show==1)
        ostruct.savename = sprintf('Fig2-IV+Ba-%s',neuron.experiment);
    else
        ostruct.savename = sprintf('Fig2-IV+Ba-%s_onlyexp',neuron.experiment);
    end
end
ostruct.handles = aGC_IVplot(expcat(targetfolder_data,'Exp_Kir',neuron_orig.experiment),targetfolder_results,ostruct);

savename = ostruct.savename ;
if ~ostruct.newborn
    ostruct.show = 2;  % 1 = show experimental data, 2 = show simulation data, 1:2 = show both
    if ~all(ostruct.show==1)
        ostruct.savename = sprintf('Fig2-IV_dyn-%s',neuron.experiment);
    else
        ostruct.savename = sprintf('Fig2-IV_dyn-%s_onlyexp',neuron.experiment);
    end
    aGC_IVplotdyn(expcat(targetfolder_data,'Exp_Kir',neuron.experiment),targetfolder_results,ostruct);
end

neuron = blockchannel(neuron,{'Kir21','pas'},[99 30]);
ostruct.bablock = 1;
aGC_IV(neuron,tree,params,targetfolder_data,ostruct);

if ostruct.newborn
    ostruct.dataset =0;
else
    ostruct.dataset =5;
end
ostruct.savename = savename;
aGC_IVplot(expcat(targetfolder_data,'Exp_Kir',neuron.experiment),targetfolder_results,ostruct)
ostruct = rmfield(ostruct,'savename');
ostruct.handles = [];



%% Mongiat FI + Ba, simulate the F-I relationship with and without blocking Kir by application of Barium, Figure 2 & 6
neuron = neuron_orig;
neuron.experiment = strcat(neuron.experiment,'_Final');
ostruct.figureheight = 4;
ostruct.figurewidth = 6;
ostruct.amp = 0:5:120;% current steps in pA which are to be simulated
ostruct.variabledt = 0;  % 1 = variable dt
ostruct.coarse = 0.5;  % 0 = dt of 0.025, 0.5 = dt of 0.05, 1 = dt of 0.1 and nseg = 1
if ostruct.newborn
    ostruct.dataset = 2.28;
    ostruct.savename = sprintf('Fig6_Final-%s',neuron.experiment);
    ostruct.savename2 = sprintf('Fig6-FI-%s_young',neuron.experiment);
    ostruct.savename3 = sprintf('Fig6-FI+Ba-%s_young',neuron.experiment);
    steps = [10, 50]/1000;  % current step [pA] of which voltage trace is plotted
else
    ostruct.dataset = 3;
    ostruct.savename = sprintf('Fig2_Final-%s',neuron.experiment);
    ostruct.savename2 = sprintf('Fig2-FI-%s',neuron.experiment);
    ostruct.savename3 = sprintf('Fig2-FI+Ba-%s',neuron.experiment);
    if isnan(ostruct.vmodel)
        steps = [30, 90]/1000;  % current step [pA] of which voltage trace is plotted
    else
        steps = [30, 75]/1000;  % current step [pA] of which voltage trace is plotted
    end
end

aGC_currsteps(neuron,tree,params,targetfolder_data,ostruct)

ostruct.show = 1:3; % 1 = show experiment, 2 = show data and 3 = make extra figures only with the current steps defined by "steps"
aGC_currsteps_plot(targetfolder_data,targetfolder_results,neuron,ostruct,steps)
% % this part is only used if experimental data and convex hulls of data
% % should be plotted (for publicatoin)
% ostruct.show = 1;
% ostruct.savename = strcat(ostruct.savename,'_onlyexp');
% [prop1,fig] =     aGC_APprop(targetfolder_data,targetfolder_results,neuron,ostruct);
% if all(ostruct.show == 1)
%     ostruct.blur = 1;
%     ostruct.patch = 1;
%     figure(fig(6))
%     [~,fig2] = plot2dens(gca,[2 2 1],ostruct);
%     figure(fig2(4))
%     tprint(fullfile2(targetfolder_results,strcat(ostruct.savename,'-','APvthreshVSAPwidth_densitypatch')),'-pdf');
%     figure(fig(7))
%     [~,fig2] = plot2dens(gca,[20 0.2 1],ostruct);
%     figure(fig2(4))
%     tprint(fullfile2(targetfolder_results,strcat(ostruct.savename,'-','APampVSfahp_densitypatch')),'-pdf');
% end
% ostruct.savename = ostruct.savename(1:end-8);

ostruct.show = 2;   % 1 = show experimental data, 2 = show simulation data, 1:2 = show both
[prop2,fig] =     aGC_APprop(targetfolder_data,targetfolder_results,neuron,ostruct,tree);

ostruct.show = 1:2;  % 1 = show experimental data, 2 = show simulation data, 1:2 = show both
ostruct.handles = [];
ostruct.savename = ostruct.savename2;
if all(ostruct.show == 1)
    ostruct.savename = strcat(ostruct.savename,'_onlyexp');
end
if ~ostruct.newborn
    ostruct.figurewidth = 4;
end
ostruct.handles = aGC_FIplot(targetfolder_data,targetfolder_results,neuron,ostruct);

if ~ostruct.newborn
    neuron = blockchannel(neuron,{'Kir21','pas'},[99 30]);
    if ostruct.newborn
        ostruct.dataset =0;
    else
        ostruct.dataset =5;
    end
    
    aGC_currsteps(neuron,tree,params,targetfolder_data,ostruct)
    
    ostruct.handles = [];
    ostruct.savename = ostruct.savename3;
    if all(ostruct.show == 1)
        ostruct.savename = strcat(ostruct.savename,'_onlyexp');
    end
    aGC_FIplot(targetfolder_data,targetfolder_results,neuron,ostruct);
    ostruct.handles = [];
    ostruct.savename = [];
end


%% Mongiat dV / dt curve,Figure 2 & 6, phase plots, extra simulation needed because a small dt is necessary to correctly assess maximal dV/dt's
neuron = neuron_orig;

neuron.experiment = strcat(neuron.experiment,'_dV');
ostruct.amp = [40,65,90,115]; % steps that are simulated plotted [pA]
ostruct.duration = 200;
ostruct.variabledt = 0;
ostruct.coarse = 0;
ostruct.show = 1:2;
ostruct.figurewidth = 6;
ostruct.figureheight = 4;
ostruct.ampprop = 90;  % plot this current step in the single figure
if ostruct.newborn
    ostruct.dataset =2.28;
    ostruct.savename = sprintf('Fig6_-%s',neuron.experiment);
else
    ostruct.dataset =3;
    ostruct.savename = sprintf('Fig2_-%s',neuron.experiment);
end
if ostruct.usemorph >= 4  % = rat morphology
    ostruct.show = 2;
    ostruct.holding_voltage = -80;
    ostruct.savename = sprintf('SupplFigX_-%s',neuron.experiment);
end

aGC_currsteps(neuron,tree,params,targetfolder_data,ostruct)

maxdv = aGC_dVplot(targetfolder_data,targetfolder_results,neuron,ostruct);


%% Kv1.1 overexpression, Figure 5, reproduce Kirchheim et al 2013, Kv1.1 overexpression after status epilepticus reduced FI and increases spiking delay
steps = [70]/1000;
neuron = neuron_orig;
ostruct.handles = [];
neuron.experiment = strcat(neuron.experiment,'_Final');
ostruct.figureheight = 4;
ostruct.figurewidth = 6;
ostruct.duration = 200;
ostruct.amp = 10:10:120; % pA
ostruct.variabledt = 0;
ostruct.coarse = 0.5;

ostruct.dataset =3;
ostruct.savename1 = sprintf('Fig5_Kv11-%s',neuron.experiment);
ostruct.savename2 = sprintf('Fig5-FI-Kv11-%s',neuron.experiment);

ostruct.show = 2:3; %& explicit spikes
ostruct.handles1 = aGC_currsteps_plot(targetfolder_data,targetfolder_results,neuron,ostruct,steps);

ostruct.show = 2;

ostruct.handles2 = aGC_FIplot(targetfolder_data,targetfolder_results,neuron,ostruct);

%increase Kv1.1 10fold
for t = 1:numel(neuron.mech)
    fields = fieldnames(neuron.mech{t});
    for f1 = 1:numel(fields)
        if isfield(neuron.mech{t}.(fields{f1}),'Kv11')
            neuron.mech{t}.(fields{f1}).Kv11.gkbar = neuron.mech{t}.(fields{f1}).Kv11.gkbar * 10;
        end
    end
end
neuron.experiment = strcat(neuron.experiment,'_Kv11');
tree_orig = tree;
for t=1:numel(tree)
    tree{t}.col = colorme('red');
end
aGC_currsteps(neuron,tree,params,targetfolder_data,ostruct)

ostruct.savename = ostruct.savename1;
ostruct.show = 2:3; %& explicit spikes
ostruct.handles = ostruct.handles1;
aGC_currsteps_plot(targetfolder_data,targetfolder_results,neuron,ostruct,steps);

ostruct.show = 2;
ostruct.savename = ostruct.savename2;

ostruct.handles = ostruct.handles2;
aGC_FIplot(targetfolder_data,targetfolder_results,neuron,ostruct);
    ostruct.handles = [];
    ostruct.savename = [];
tree = tree_orig;

%% Synaptic stimulation, Figure 7
ostruct.handles = [];
type = 'temporal'; % test, theta-burst stimulation, regular input, temporal shift in input, spatial shift in input
ostruct.synmode = 1; % 1 = no adjustment of synapse number, 2 = newborns have only half of synapse number, 3 = adjust the number of synapses to get subthreshold response
% ostruct.figureheight = 3;
aGC_synstim(params,neuron,tree,type,ostruct,targetfolder_data)

ostruct.handles = aGC_synstim_plot(params,tree,type,ostruct,targetfolder_data,targetfolder_results,0);

%% *********************************  RAT EXPERIMENTS **********************************

%% passive parameter assessment rat
ostruct.passtest = 'Std'; % different protocols to measure Rin/tau/capacitance/Vrest (not all can do everything): Mongiat Mongiat2 SH Brenner Std
[Rin, tau, cap, Vrest] = aGC_passtests(neuron,tree,params,targetfolder_results,ostruct);
%% IV rat, generate I-V relationship curve and additinally plot Rin measurements from rat, Figure 3
ostruct.holding_voltage = -80; % mV
ostruct.subtract_hv = 1; % boolean subtract holding voltage current
ostruct.show = 2; % 0= nothing, 1 = only exp data, 2 = only model data
ostruct.variabledt = 1; % 1 = variable dt
ostruct.coarse = 0;
ostruct.single = 0; % show single data curves instead of mean

aGC_IV(neuron,tree,params,targetfolder_data,ostruct);

ostruct.savename = sprintf('Fig3-IV_Mehranfard15-%s',neuron.experiment);
if ostruct.ratadjust
    ostruct.savename = strcat(ostruct.savename,'_ratadjust');
end

aGC_IVplot(expcat(targetfolder_data,'Exp_Kir',neuron.experiment),targetfolder_results,ostruct);
% aGC_IVplotdyn(expcat(targetfolder_data,'Exp_Kir',neuron.experiment),targetfolder_results,ostruct);


%% Rat FI curve, Figure 3, similar to Mehranfard 2015 Plos One 
neuron = neuron_orig;
if isfield(ostruct,'handles')
    ostruct = rmfield(ostruct,'handles');
end
ostruct. handles = [];
% params.celsius = 24;
ostruct.amp = 50:50:300; % pA
ostruct.duration = 1000;
ostruct.variabledt = 0;
ostruct.coarse = 0.5;
ostruct.holding_voltage = -80;  % leider unknown. ABER OHNE LJP
% ostruct = rmfield(ostruct,{'handles'});

aGC_currsteps(neuron,tree,params,targetfolder_data,ostruct)

% ostruct.savename = sprintf('Fig3_spiking_Mehranfahrd15-%s',neuron.experiment);
% if ostruct.ratadjust
%     ostruct.savename = strcat(ostruct.savename,'_ratadjust');
% end
ostruct.dataset = 0;  % 1 = old mature dataset, 2 = old young dataset, 3 = new mature dataset, 4 = new young dataset, 5 = new mature BaCl dataset, 6 = new young BaCl dataset
ostruct.show = 1:2; % data (Meranfahrd) &  simulation
ostruct.savename = sprintf('Fig3_FI_Mehranfahrd15-%s',neuron.experiment);
if ostruct.ratadjust
    ostruct.savename = strcat(ostruct.savename,'_ratadjust');
end
aGC_FIplot(targetfolder_data,targetfolder_results,neuron,ostruct);
ostruct.show = 2:3; %   simulation & explicit spikes
ostruct.savename = sprintf('Fig3_Final-%s',neuron.experiment);
if ostruct.ratadjust
    ostruct.savename = strcat(ostruct.savename,'_ratadjust');
end
aGC_currsteps_plot(targetfolder_data,targetfolder_results,neuron,ostruct,0.1);

%% ISI spike adaptaption Figure 5, similar to MA14 Fig 8
%!!!
params.celsius = 33;
%!!!
neuron = neuron_orig;
ostruct.holding_voltage = -62;
ostruct.find_freq = 6; % 6 spikes
neuron.experiment = strcat(sprintf('MA14Fig8_freq%g_%dC_',ostruct.find_freq,params.celsius),neuron.experiment);

if isfield(ostruct,'handles')
    ostruct = rmfield(ostruct,'handles');
end
ostruct. handles = [];
ostruct.amp = 200;%pA ACHTUNG, das wird an sich ignoriert durch das find_freq!
steps = [0.15];%,0.2];
ostruct.coarse = 0.5;
ostruct.show = 2; % only simulation
ostruct.duration = 100;
ostruct.variabledt = 0;
ostruct.ampprop = 200; 
ostruct.savename = 'Fig5-ISIadapt100ms_MA14';
if ostruct.ratadjust
    ostruct.savename = strcat(ostruct.savename,'_ratadjust');
end
ostruct.dataset = 0;  % 1 = old mature dataset, 2 = old young dataset, 3 = new mature dataset, 4 = new young dataset, 5 = new mature BaCl dataset, 6 = new young BaCl dataset

aGC_currsteps(neuron,tree,params,targetfolder_data,ostruct)

aGC_currsteps_plot(targetfolder_data,targetfolder_results,neuron,ostruct,steps);

ostruct.show = 0;
props = aGC_APprop(targetfolder_data,targetfolder_results,neuron,ostruct);

ISIadp{1} = NaN(numel(tree),1);
for t = 1:numel(tree)
    if numel(props.APISI{1,t}) > 1
        ISIadp{1}(t) = 1-props.APISI{1,t}(1)/props.APISI{1,t}(end);
    end
end
ostruct.show = 2;

neuron = blockchannel(neuron,{'SK2'},[100]);
% ostruct.dataset =5;
aGC_currsteps(neuron,tree,params,targetfolder_data,ostruct)
aGC_currsteps_plot(targetfolder_data,targetfolder_results,neuron,ostruct,steps);

ostruct.show = 0;
props = aGC_APprop(targetfolder_data,targetfolder_results,neuron,ostruct);

ISIadp{2} = NaN(numel(tree),1);
for t = 1:numel(tree)
    if numel(props.APISI{1,t}) > 1
        ISIadp{2}(t) = 1-props.APISI{1,t}(1)/props.APISI{1,t}(end);
    end
end

ostruct.show = 2;

neuron = neuron_orig;
neuron.experiment = strcat(sprintf('MA14Fig8_freq%g_%dC_',ostruct.find_freq,params.celsius),neuron.experiment);
neuron = blockchannel(neuron,{'Kv723'},[100]);
% ostruct.dataset =5;
aGC_currsteps(neuron,tree,params,targetfolder_data,ostruct)
aGC_currsteps_plot(targetfolder_data,targetfolder_results,neuron,ostruct,steps);

ostruct.show = 0;
props = aGC_APprop(targetfolder_data,targetfolder_results,neuron,ostruct);

ISIadp{3} = NaN(numel(tree),1);
for t = 1:numel(tree)
    if numel(props.APISI{1,t}) > 1
        ISIadp{3}(t) = 1-props.APISI{1,t}(1)/props.APISI{1,t}(end);
    end
end
ostruct.show = 2;

dstruct = [];gstruct = [];
dstruct.data.ISIadapt = [];
dstruct.strs.ISIadapt.ylab = 'Adaptation index [%]';
gstruct.group{1} = [ones(1,numel(tree)), repmat(2,1,numel(tree)), repmat(3,1,numel(tree))];
gstruct.ugroupdef{1} = {'CTRL','SK block','Kv7 block'};
dstruct.strs.ISIadapt.xstr = gstruct.ugroupdef{1};

for n = 1:3
    dstruct.data.ISIadapt(1,(n-1)*numel(tree)+1:numel(tree)*n) = ISIadp{n}./ISIadp{1}*100; 
end

handle = Megaplotter(dstruct,gstruct,[],ostruct);
neuron = neuron_orig;
neuron.experiment = strcat(sprintf('MA14Fig8_freq%g_%dC_',ostruct.find_freq,params.celsius),neuron.experiment);
tprint(fullfile2(targetfolder_results,expcat('Fig.5-ISIadapt',neuron.experiment)),'-pdf')

params.celsius = 24;
ostruct = rmfield(ostruct,'find_freq');
%% bAP simulation (Krueppel 2011) +  Calcium dynamics (Stocca 2008) Teil Ratte, Figure 4
neuron = neuron_orig;
ostruct.simple = 0;
ostruct.reduce = 0;
ostruct.dist = 'Eucl.'; % PL., Eucl.
ostruct.relamp = 0;
celsius_orig = params.celsius;
params.celsius = 33;
neuron = Q10pas(neuron,params.celsius);
neuron.experiment = sprintf('%s_%d°',neuron.experiment,params.celsius);
ostruct.show = 1;

aGC_bAP(neuron,tree,params,targetfolder_data,ostruct)


params.celsius = celsius_orig; % back to normal temp
aGC_bAP_plot(targetfolder_data,targetfolder_results,neuron,ostruct);

%% AP width BK/Kv34 + spike adaptation, Figure 5
neuron = neuron_orig;
neuron.experiment = strcat(neuron.experiment,'_APwidth');
% neuron.experiment = strcat('APwidth+adapt_',neuron.experiment);
ostruct.amp = [90,250]; % pA  < 10 HZ and höher
ostruct.duration = 1000;
ostruct.variabledt = 0;
ostruct.coarse = 0.5;
ostruct.holding_voltage = -80;
ostruct.ampprop = 90;
% ostruct = rmfield(ostruct,{'handles'});
if ostruct.newborn
    ostruct.savename = 'Fig5-APprop_newborn';
else
    ostruct.savename = 'Fig5-APprop';
end

aGC_currsteps(neuron,tree,params,targetfolder_data,ostruct)

ostruct.handles = [];
ostruct.show = 0;

prop = [];
for a = 1:2
    ostruct.ampprop = ostruct.amp(a);
    prop{a}(1) = aGC_APprop(targetfolder_data,targetfolder_results,neuron,ostruct);
end


ostruct.show = 2;

ostruct.handles = aGC_currsteps_plot(targetfolder_data,targetfolder_results,neuron,ostruct);
for t = numel(tree):-1:1
    if t == 1
        ostruct.handles(1).Children(end-1).Children(t).Color = [0 0 0];
        ostruct.handles(1).Children(end).Children(t).Color = [0 0 0];
    else
        delete(ostruct.handles(1).Children(end-1).Children(t))
        delete(ostruct.handles(1).Children(end).Children(t))
    end
end

neuron = blockchannel(neuron,{'Kv34'},100);
%
aGC_currsteps(neuron,tree,params,targetfolder_data,ostruct)
%
ostruct.show = 0;
for a = 1:2
    ostruct.ampprop = ostruct.amp(a);
    prop{a}(2) = aGC_APprop(targetfolder_data,targetfolder_results,neuron,ostruct);
end
ostruct.show = 2;
aGC_currsteps_plot(targetfolder_data,targetfolder_results,neuron,ostruct);
for t = numel(tree):-1:1
    if t == 1
        ostruct.handles(1).Children(end-1).Children(t).Color = [1 0 0];
        ostruct.handles(1).Children(end).Children(t).Color = [1 0 0];
    else
        delete(ostruct.handles(1).Children(end-1).Children(t))
        delete(ostruct.handles(1).Children(end).Children(t))
    end
end

neuron = neuron_orig;
neuron.experiment = strcat(neuron.experiment,'_APwidth');
neuron = blockchannel(neuron,{'BK'},100);
%
aGC_currsteps(neuron,tree,params,targetfolder_data,ostruct)
%
ostruct.show = 0;
for a = 1:2
    ostruct.ampprop = ostruct.amp(a);
    prop{a}(3) = aGC_APprop(targetfolder_data,targetfolder_results,neuron,ostruct);
end
ostruct.show = 2;
aGC_currsteps_plot(targetfolder_data,targetfolder_results,neuron,ostruct);
for t = numel(tree):-1:1
    if t == 1
        ostruct.handles(1).Children(end-1).Children(t).Color = [0 0 1];
        ostruct.handles(1).Children(end).Children(t).Color = [0 0 1];
    else
        delete(ostruct.handles(1).Children(end-1).Children(t))
        delete(ostruct.handles(1).Children(end).Children(t))
    end
end
% % SK/Kv7 block
% 
% %
dstruct = [];gstruct = [];
dstruct.data.APwidth = [];
dstruct.data.APwidth2 = [];
dstruct.strs.APwidth.xstr = sprintfc('%d pA',ostruct.amp );
dstruct.strs.APwidth.ylab = 'Half-amplitude width [ms]';
dstruct.strs.APwidth2.xstr = sprintfc('%d pA',ostruct.amp );
dstruct.strs.APwidth2.ylab = 'Spike width [ms]';
gstruct.group{1} = [ones(1,numel(tree)), repmat(2,1,numel(tree)), repmat(3,1,numel(tree))];%[4,4,5,5,6,6,ones(1,numel(tree)), repmat(2,1,numel(tree)), repmat(3,1,numel(tree))];
gstruct.ugroupdef{1} = {'CTRL','Kv3.4 block','BK block'};%,'CTRL data','Kv3.4 block data','BK block data'};
xth = [1,5];
for a = 1:2
%     ostruct.handles(1).Children(end-2+a).XLabel.String = 'Time [ms]';
%     ostruct.handles(1).Children(end-2+a).YLabel.String = 'Membrane Potential [mV]';
    for n = 1:3
        dstruct.data.APwidth(a,(n-1)*numel(tree)+1:numel(tree)*n) = cellfun(@(x) x(xth(a)),prop{a}(n).APwidth); % get AP widths of xth AP in all simulations
        dstruct.data.APwidth2(a,(n-1)*numel(tree)+1:numel(tree)*n) = cellfun(@(x) x(xth(a)),prop{a}(n).APwidth2); % get AP widths of xth AP in all simulations
    end
end
additionalData = [0.9662-0.0185,0.9662+0.0185,NaN,NaN,1.5748-0.0154,1.5748+0.0154;1.5345+0.1081,1.5345-0.1081,NaN,NaN,2.1089-0.3732,2.1089+0.3732];
additionalData(:,3:4) = [mean(additionalData(:,1:2),2)*(1.0793-0.0205),mean(additionalData(:,1:2),2)*(1.0793+0.0205)];
additionalData = cat(2,additionalData(:,1:2),NaN(2,numel(tree)-2),additionalData(:,3:4),NaN(2,numel(tree)-2),additionalData(:,5:6),NaN(2,numel(tree)-2));
dstruct.data.APwidth = cat(1,additionalData,dstruct.data.APwidth);
dstruct.data.APwidth2 = cat(1,additionalData,dstruct.data.APwidth2);


ostruct.handles(1).Children(end).Children(3).XData = ostruct.handles(1).Children(end).Children(3).XData - prop{1}(1).APt{end}(xth(1));
ostruct.handles(1).Children(end).Children(2).XData = ostruct.handles(1).Children(end).Children(2).XData - prop{1}(2).APt{end}(xth(1));
ostruct.handles(1).Children(end).Children(1).XData = ostruct.handles(1).Children(end).Children(1).XData - prop{1}(3).APt{end}(xth(1));

ostruct.handles(1).Children(end-1).Children(3).XData = ostruct.handles(1).Children(end-1).Children(3).XData - prop{2}(1).APt{end}(xth(2));
ostruct.handles(1).Children(end-1).Children(2).XData = ostruct.handles(1).Children(end-1).Children(2).XData - prop{2}(2).APt{end}(xth(2));
ostruct.handles(1).Children(end-1).Children(1).XData = ostruct.handles(1).Children(end-1).Children(1).XData - prop{2}(3).APt{end}(xth(2));

ostruct.handles(1).Children(end).YLim = [-80 70];
ostruct.handles(1).Children(end).XLim = [-5 20];
ostruct.handles(1).Children(end-1).YLim = [-80 70];
ostruct.handles(1).Children(end-1).XLim = [-5 20];
FontResizer
FigureResizer(5,8)
% if ostruct.newborn
%     tprint(fullfile2(targetfolder_results,'RobustnessMatrix_newborn'),'-pdf');
% else
    tprint(fullfile2(targetfolder_results,sprintf('Fig.5_APwidthSpiking_%s',neuron.experiment)),'-pdf');
    
% end
ostruct.gap = 0.3;
handle = Megaplotter(dstruct,gstruct,[],ostruct);
figure(handle{1})
tprint(fullfile2(targetfolder_results,sprintf('Fig.5_APhalfwidth_%s',neuron.experiment)),'-pdf');
figure(handle{2})
tprint(fullfile2(targetfolder_results,sprintf('Fig.5_APwidth_%s',neuron.experiment)),'-pdf');

%% resonance test with Ba and ZD application, Suppl. Fig.

neuron = neuron_orig;
neuron.experiment = strcat(neuron.experiment,'_resonance');
aGC_resonance(params,neuron,tree,ostruct,targetfolder_results)


%% Sensitivity Matrix Fig. 5 & Suppl. Fig. 
type = 'up'; % possible values: up down   ... regulate channels up or down
neuron = neuron_orig;
params = params_orig;
params.exchfolder = 't2nexchange_aGCmorphsim4';
neuron.experiment = strcat('robust_',neuron.experiment);
changs = {'','cm','Ra','pas','Kir21','HCN','cAMP','na8st','Kv14','Kv21','Kv34','Kv42','Kv7','BK','SK','Cav12','Cav13','Cav22','Cav32','E_K','E_N_a','E_P_a_s','[Ca]o','temp'};
% rmatrix = zeros(15,numel(changs));
ostruct.show = 0;

if strcmp(type,'up')
    amount = [-100,-10,20,4];
else
    amount = [50,10,-20,1]; 
end
w = waitbar(0,'Matrix is calculating, please get yourself a snickers...');
for v = 1:numel(changs) %
    
    switch changs{v}
%         case {'pas','na8st','Kir21','HCN','Kv14','Kv21','Kv34','Kv42','Kv7','SK','Cav12','Cav13','Cav32','BK'}
            
%         case {'gabkbar'}
%             neuron = blockchannel(neuron,'BK',50,changs{v}); % specify BK channel
        case 'cm'
            neuron = blockchannel(neuron,'pas',amount(1),changs{v}); % specify pas channel
        case 'Ra'
            neuron = blockchannel(neuron,'pas',amount(1),changs{v}); % specify pas channel
        case 'E_K'
            for t = 1:numel(tree)
                neuron.mech{t}.all.k_ion.ek = neuron.mech{t}.all.k_ion.ek + amount(2);
            end
            neuron.experiment = strcat(neuron.experiment,'_ek');
        case 'E_N_a'
            for t = 1:numel(tree)
                neuron.mech{t}.all.na_ion.ena = neuron.mech{t}.all.na_ion.ena + amount(3);
            end
            neuron.experiment = strcat(neuron.experiment,'_ena');
        case 'E_P_a_s'
            for t = 1:numel(tree)
                fields = fieldnames(neuron.mech{t});
                for f = 1:numel(fields)
                    if isfield(neuron.mech{t}.(fields{f}),'pas')
                        neuron.mech{t}.(fields{f}).pas.e = neuron.mech{t}.(fields{f}).pas.e + amount(2);
                    end
                end
            end
            neuron.experiment = strcat(neuron.experiment,'_epas');
        case '[Ca]o'
            for t = 1:numel(tree)
                neuron.mech{t}.all.ca_ion.cao0 = amount(4);%neuron.mech{t}.all.ca_ion.cao0 / 2;
            end
            neuron.experiment = strcat(neuron.experiment,'_eca');
        case 'temp'
            params.celsius = params.celsius + amount(2);
            neuron.experiment = strcat(neuron.experiment,'_temp');
        case 'cAMP'
                neuron = changecAMP(neuron,1); % 1µM cAMP
        otherwise
            if ~isempty(changs{v})
                
                neuron = blockchannel(neuron,changs{v},amount(1));
            end
%             errordlg('not found')
%             return
           
    end
    %a
    ostruct.passtest = 'Std'; %Mongiat Mongiat2 SH Brenner Std
    [Rin, tau, ~, Vrest] = aGC_passtests(neuron,tree,params,targetfolder_results,ostruct);
    % caution, VoltageClamps add the series resistance of the electrode (10MOhm
    % by now). With passive, capacitance is higher..? But Rin is always same
    rmatrix(1,v) = mean(Vrest);
    rmatrix(2,v) = mean(Rin);
    rmatrix(3,v) = mean(tau);
    % b
    ostruct.holding_voltage = -80; %  % CAUTION! no LJP subtraction if concrete value is given
    ostruct.duration = 200;%200;%200;
    ostruct.amp = 90; %15:5:90 pA  brauch die steps NICHT MEHR
    ostruct.variabledt = 0;
    ostruct.coarse = 0.5;
    ostruct.ampprop = 90;
    ostruct.data = 2;
    rmatrix(5,v) = nanmean(find_curr(params,neuron,tree,'spike'));
    
    aGC_currsteps(neuron,tree,params,targetfolder_data,ostruct)
    props = aGC_APprop(targetfolder_data,targetfolder_results,neuron,ostruct);
    rmatrix(4,v) = nanmean(cellfun(@(y) nanmean(y),props.APiv(1,:)));
%     rmatrix(5,v) = nanmean(cellfun(@(y) nanmean(y),props.APic(1,:)));
    rmatrix(6,v) = nanmean(cellfun(@(y) nanmean(y),props.APamp(1,:)));
    rmatrix(7,v) = nanmean(cellfun(@(y) nanmean(y),props.APwidth(1,:)));
    rmatrix(8,v) = nanmean(cellfun(@(y) nanmean(y),props.fAHPabs(1,:)));
    % c
    ostruct.duration = 1000;%200;%200;
    ostruct.amp = [100,200]; % pA
    ostruct.ampprop = 100;
    aGC_currsteps(neuron,tree,params,targetfolder_data,ostruct)
    props = aGC_APprop(targetfolder_data,targetfolder_results,neuron,ostruct);
    rmatrix(9,v) =  mean(cellfun(@(y) numel(y),props.APind(1,:))); %# APs@100
    rmatrix(11,v) =  nanmean(cellfun(@(y) mean(y),props.APISI(1,:))); % ISI mean
%     try
%         rmatrix(12,v) =  mean(cellfun(@(y) 1-(y(1)/y(end)),props.APISI(1,:)));% ISI adaptation
%     catch   % possibly no spike
       isia = NaN(numel(tree),1);
       for t = 1:numel(tree)
           if numel(props.APISI{1,t}) > 1
               isia(t) = 1-props.APISI{1,t}(1)/props.APISI{1,t}(end);
               if isia(t) == 0
                   isia(t) = NaN;
               end
           end
       end
       rmatrix(12,v) = nanmean(isia);
%     end ISIadp{3}(t) = 1-props.APISI{1,t}(1)/props.APISI{1,t}(end);
%     aGC_currsteps_plot(targetfolder_data,targetfolder_results,neuron,ostruct);

    %d
    ostruct.ampprop = 200;
    props = aGC_APprop(targetfolder_data,targetfolder_results,neuron,ostruct);
    rmatrix(10,v) = mean(cellfun(@(y) numel(y),props.APind(1,:)));%# APs@200
    %e
    ostruct.simple = 0;
    ostruct.reduce = 1;
    ostruct.relamp = 0;
    aGC_bAP(neuron,tree,params,targetfolder_data,ostruct)
    [bAPdisthm,mveloc_dend,mveloc_farax,mveloc_nearax] = aGC_bAP_plot(targetfolder_data,targetfolder_results,neuron,ostruct);
    rmatrix(13,v) = nanmean(mveloc_dend);
    rmatrix(14,v) = nanmean(mveloc_farax);
%     rmatrix(15,v) = nanmean(bAPdisthm);
    
    neuron = neuron_orig;
    params = params_orig;
    neuron.experiment = strcat('robust_',neuron.experiment);
    params.exchfolder = 't2nexchange_aGCmorphsim4';
    waitbar(v/numel(changs),w)
end
close(w)
values = {'V_R_e_s_t','R_i_n','membr. time const.','voltage thresh. AP','current thresh. AP','AP amp','AP width','AP fAHP','#APs  100 pA','#APs  200 pA','ISI','adaptation ratio','dendr. veloc  mean','ax. veloc  far ax'};%,'bAP PL dist half-max'};

figure;
if any(imag(rmatrix(:)))
   [x,y] = find(imag(rmatrix));
   errordlg(sprintf('Warning, found imaginary values for %s during change of %s',values{x},changs{y}))
end
imagesc(real(rmatrix(:,2:end))./repmat(real(rmatrix(:,1)),1,size(rmatrix,2)-1)-1)
set(gca,'YTick',1:15)
set(gca,'YTickLabel',values)
set(gca,'XTick',1:numel(changs)-1)
set(gca,'XTickLabel',changs(2:end))

set(gca,'CLim',[-0.5 0.5])
set(gca,'Box','off')
colorbar
o.image = 0;
o.border = [8,2];
FigureResizer(8,20,[],o)
if ostruct.newborn
    save(fullfile(targetfolder_data,sprintf('Fig.5_SensitivityMatrix_%s_newborn_%s.mat',type,neuron.experiment)),'rmatrix')
    tprint(fullfile2(targetfolder_results,sprintf('Fig.5_SensitivityMatrix_%s_newborn_%s',type,neuron.experiment)),'-pdf');
else
    save(fullfile(targetfolder_data,sprintf('Fig.5_SensitivityMatrix_%s_%s.mat',type,neuron.experiment)),'rmatrix')
    tprint(fullfile2(targetfolder_results,sprintf('Fig.5_SensitivityMatrix_%s_%s',type,neuron.experiment)),'-pdf');
end


%% ********** NOT USED IN PUBLICATION ****************
%% Brenner FI normal, this is the strong mouse stimulation protocol of Brenner et al 2005
neuron = neuron_orig;
ostruct.amp = 50:50:400; % pA
ostruct.duration = 900;
ostruct.variabledt = 0;
ostruct.coarse = 0.5;
ostruct.holding_voltage = -80;
ostruct.show = 1:2;
if ostruct.newborn
    ostruct.savename = sprintf('Fig2_FI_Brenner_newborn-%s',neuron.experiment);
else
    ostruct.savename = sprintf('Fig2_FI_Brenner-%s',neuron.experiment);
end
ostruct.dataset = 7;  % 1 = old mature dataset, 2 = old young dataset, 3 = new mature dataset, 4 = new young dataset, 5 = new mature BaCl dataset, 6 = new young BaCl dataset, 7 = Brenner mouse
aGC_currsteps(neuron,tree,params,targetfolder_data,ostruct)

aGC_FIplot(targetfolder_data,targetfolder_results,neuron,ostruct);
aGC_currsteps_plot(targetfolder_data,targetfolder_results,neuron,ostruct)

%% attempt to put together real and synth results for bAP
nneuron{1}.experiment = 'rat_mGC_Beining_scaledspines_ratadjust_Q10pas_33°';
nneuron{2}.experiment = 'rat_matGC_art_scaledspines_ratadjust_Q10pas_33°';
params.celsius = celsius_orig; % back to normal temp
aGC_bAP_plot(targetfolder_data,targetfolder_results,nneuron,ostruct);
%% Calcium channel contributions (Eliot Johnston 1994), no Figure
neuron = neuron_orig;

aGC_Ca_proportions_Elliot(neuron,tree,params,targetfolder_results)

%% some FI curves similar to Mateos-Aparicio 2014
%!!!
params.celsius = 33;
%!!!
neuron = neuron_orig;
neuron.experiment = strcat(sprintf('MA14_%dC_',params.celsius),neuron.experiment);

if isfield(ostruct,'handles')
    ostruct = rmfield(ostruct,'handles');
end
ostruct. handles = [];
ostruct.amp = 50:50:200;%50:50:350; % pA
steps = [0.1,0.2];
ostruct.coarse = 0.5;
ostruct.show = 2; % only simulation
ostruct.duration = 1000;
ostruct.variabledt = 0;
ostruct.holding_voltage = -77;
ostruct.ampprop = 100; % necessary for instFI
% ostruct = rmfield(ostruct,{'handles'});
ostruct.savename = 'Fig5-instFI_MA14';
if ostruct.ratadjust
    ostruct.savename = strcat(ostruct.savename,'_ratadjust');
end
ostruct.dataset = 0;  % 1 = old mature dataset, 2 = old young dataset, 3 = new mature dataset, 4 = new young dataset, 5 = new mature BaCl dataset, 6 = new young BaCl dataset

aGC_currsteps(neuron,tree,params,targetfolder_data,ostruct)
% 
% % ostruct.handles = [];
% ostruct.handles = aGC_FIplot(targetfolder_data,targetfolder_results,neuron,ostruct);
% % ostruct.handles = 
aGC_currsteps_plot(targetfolder_data,targetfolder_results,neuron,ostruct,steps);
% ostruct.handles(1).Children(end).Children(1).Color = [0 0 0];
% ostruct.handles(2).Children(end).Children(1).Color = [0 0 0];
% ostruct.handles(3).Children(end).Children(1).Color = [0 0 0];
% figure(ostruct.handles(2))
% dataMA = importdata(fullfile2(targetfolder_results,'FI_MA-CTRL.csv')); 
% plot (dataMA(:,1)*1000, dataMA(:,2), 'k')
% 
% tmphandles = ostruct.handles; ostruct.handles =[];
% ostruct.handles2 = aGC_currsteps_plot(targetfolder_data,targetfolder_results,neuron,ostruct);
% ostruct.handles = tmphandles;

neuron = blockchannel(neuron,{'SK2'},[100]);
% ostruct.dataset =5;
aGC_currsteps(neuron,tree,params,targetfolder_data,ostruct)
aGC_currsteps_plot(targetfolder_data,targetfolder_results,neuron,ostruct,steps);
% aGC_FIplot(targetfolder_data,targetfolder_results,neuron,ostruct);
% ostruct.handles(1).Children(end).Children(1).Color = [1 0 0];
% ostruct.handles(2).Children(end).Children(1).Color = [1 0 0];
% ostruct.handles(3).Children(end).Children(1).Color = [1 0 0];
% figure(ostruct.handles(2))
% dataMA = importdata(fullfile2(targetfolder_results,'FI_MA-Apamin.csv'));
% plot (dataMA(:,1)*1000, dataMA(:,2), 'r')

% tmphandles = ostruct.handles; ostruct.handles = ostruct.handles2;
% ostruct.handles2 = aGC_currsteps_plot(targetfolder_data,targetfolder_results,neuron,ostruct);
% ostruct.handles = tmphandles;

neuron = neuron_orig;
neuron.experiment = strcat(sprintf('MA14_%dC_',params.celsius),neuron.experiment);
neuron = blockchannel(neuron,{'Kv723'},[100]);
% ostruct.dataset =5;
aGC_currsteps(neuron,tree,params,targetfolder_data,ostruct)
aGC_currsteps_plot(targetfolder_data,targetfolder_results,neuron,ostruct,steps);

% aGC_FIplot(targetfolder_data,targetfolder_results,neuron,ostruct);
% ostruct.handles(1).Children(end).Children(1).Color = [0 0.5 0];
% ostruct.handles(2).Children(end).Children(1).Color = [0 0.5 0];
% ostruct.handles(3).Children(end).Children(1).Color = [0 0.5 0];
% figure(ostruct.handles(2))
% dataMA = importdata(fullfile2(targetfolder_results,'FI_MA-XE991.csv'));
% plot (dataMA(:,1)*1000, dataMA(:,2), 'Color',[0 0.5 0])
% tprint(fullfile2(targetfolder_results,expcat('Fig.5-FI',neuron.experiment)),'-pdf')
% figure(ostruct.handles(1))
% tprint(fullfile2(targetfolder_results,expcat('Fig.5-instFI',neuron.experiment)),'-pdf')

% tmphandles = ostruct.handles; ostruct.handles = ostruct.handles2;
% ostruct.handles2 = aGC_currsteps_plot(targetfolder_data,targetfolder_results,neuron,ostruct);
% ostruct.handles = tmphandles;



params.celsius = 24;
