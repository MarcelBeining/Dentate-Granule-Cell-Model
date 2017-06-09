function aGC_sAHPstimMA14(neuron,tree,params,targetfolder)

cstep = 0.2; %nA !
params.accuracy = 1;  % for more nseg in axon and soma!
params.dt=0.05;
params.cvode = 1;
params.tstop = 350;    
params.skiprun = 0; %!!!!!!!!!
if ~exist('hstep','var')
    hstep = [];
end
hstep = t2n_findCurr(tree,params,neuron,-62,hstep,'-q-d');

for t = 1:numel(tree)
    nodes{t} = 1;
    neuron.record{t}.cell = struct('node',nodes{t},'record','v');
    neuron.pp{t}.IClamp = struct('node',1,'times',[-400,50,200],'amp', [hstep(t) hstep(t)+cstep hstep(t)]); %n,del,dur,amp
end
nneuron = cell(1,2);
for s = 1:2
    if s == 2
        nneuron{s} = t2n_t2n_blockchannel(neuron,'SK');
    else
        nneuron{s} = neuron;
    end
end
[out, minterf] = t2n(tree,params,nneuron,'-q-d-w');
% out = out{1};
if isfield(out,'error')
    return
end
save(fullfile(targetfolder,sprintf('Exp_MA14stimsAHP%s.mat',neuron.experiment)),'tree','out')