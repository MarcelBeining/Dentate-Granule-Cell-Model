function aGC_spikingadaptation(neuron,tree,params,targetfolder_data,holding_voltage)
% 
current = 150*0.001;


params.celsius = 33;
params.accuracy = 1;  % for more nseg in axon and soma!
params.v_init = -80;%-LJP;
params.tstop = 1500;
params.dt=0.05;
params.cvode = 1;
    
hstep = find_curr(params,neuron,tree,holding_voltage,[],'-q-d');

for t=1:numel(tree)
    neuron.APCount{t} = [1,-30];
end

nneuron{1} = neuron;
for t = 1:numel(tree)
    nneuron{1}.pp{t}.IClamp = struct('node',1,'times',[-100,50,1050],'amp', [hstep(t) hstep(t)+current hstep(t)]); %n,del,dur,amp
    nneuron{1}.record{t}.cell = struct('node',1,'record','v');
end


% params.skiprun = 1 %!!!!!!!!!
    [out, ~] = t2n(tree,params,nneuron,'-q-d-w');
    if out{1}.error
        return
    end

    numspikes = zeros(numel(tree),1);
    for t = 1:numel(tree)
        if isfield(out{1},'error') && out{1}.error > 0
            vol_new_curr_dend{t,s} = [] ;
            tvol_new_curr_dend{t,s} = [];
            timespikes{t} = [];
        else
            vol_new_curr_dend{t,1} = out{1}.record{t}.cell.v{1} ;
            tvol_new_curr_dend{t,1} = out{1}.t;
%             numspikes(t) = numel(out{1}.APCtimes{t}{1});
            timespikes{t} = out{1}.APCtimes{t}{1};
        end
    end
% if params.cvode
    save(sprintf(fullfile2(targetfolder_data,'Exp_Adaptation_%s.mat'),neuron.experiment),'vol_new_curr_dend','tvol_new_curr_dend','timespikes','params','current','tree','neuron')
% else
%     save(sprintf(fullfile2(targetfolder_data,'EphysModel/Exp_Spiking%s_fixed-dt.mat'),neuron.experiment),'vol_new_curr_dend','tvol_new_curr_dend','timespikes','params','cstepsSpikingModel','cstepsSpiking','tree')
end