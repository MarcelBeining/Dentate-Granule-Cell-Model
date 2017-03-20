function aGC_IV(neuron,tree,params,targetfolder_data,ostruct)
if nargin < 5 || ~isfield(ostruct,'holding_voltage')
    ostruct.holding_voltage = -80;
end
elecnode = 1;
if ~isfield(ostruct,'extract_kir')
    ostruct.extract_kir = 0;
end

LJP = params.LJP;  % use LJP only for current clamps now...
if ostruct.coarse
    vstepsModel =  (-130:10:-40) - LJP;
else
    vstepsModel =  (-130:5:-40) - LJP;
end
holding_voltage = ostruct.holding_voltage - LJP;
params.prerun = 300;
params.skiprun = 0;
dur = [105 100 105];
params.tstop = 310;
if ostruct.variabledt
    params.cvode = 1;
else
    params.dt=0.25;
    params.cvode = 0;
end

for f = 1:ostruct.extract_kir+1
    if f == 2
        neuron2 = neuron; % backup neuron struct
        neuron = blockchannel(neuron,'Kir21',100); % block Kir current
    end
    nneuron = cell(numel(vstepsModel),1);
    for s = 1:numel(vstepsModel)
        if s == 1
            nneuron{s} = neuron;
        end
        amp = cat(2,holding_voltage, vstepsModel(s), holding_voltage);       % as in Mongiat2009 (VClamp at holding potential of -70mV)
        for t = 1:numel(tree)
            nneuron{s}.pp{t}.SEClamp = struct('node',elecnode,'rs',15,'dur', dur,'amp', amp);%,'rate',1); %n,del,dur,amp
            %         nneuron{s}.record{t}.cell = struct('record','v','node',elecnode);
            nneuron{s}.record{t}.SEClamp = struct('record','i','node',elecnode);
        end
        nneuron{s} = t2n_as(1,nneuron{s});
    end
    out = t2n(tree,params,nneuron,'-q-d-w');
    if any(cellfun(@(x) x.error,out(cellfun(@(x) isfield(x,'error'),out))))
        return
    end
    if f == 1
        out2 = out;
    else
        neuron = neuron2;  % reput old neuron structure
        for s = 1:numel(vstepsModel)
            for t = 1:numel(tree)
                out{s}.record{t}.SEClamp.i{elecnode} = (interp1(out2{s}.t,out2{s}.record{t}.SEClamp.i{elecnode},0:0.5:params.tstop)-interp1(out{s}.t,out{s}.record{t}.SEClamp.i{elecnode},0:0.5:params.tstop))'; % subtract current from current with Kir blocked gives Kir current
            end
            out{s}.t = (0:0.5:params.tstop)'; % since cvode forced me to do interpolation, update t vec
        end
    end
end

ind = vstepsModel == holding_voltage;
for t = 1:numel(tree)
    mholding_current(t) = mean(out{ind}.record{t}.SEClamp.i{1}(find(out{ind}.t>=180,1,'first'):find(out{ind}.t>=200,1,'first')) *1000); % + 6 * mean(outleaksub{s}.record{t}.i{1}(find(outleaksub{s}.t>=180,1,'first'):find(outleaksub{s}.t>=200,1,'first')) *1000);
end
for s = 1:numel(vstepsModel)
    for t = 1:numel(tree)
        newcurr_dend(s,t) =  mean(out{s}.record{t}.SEClamp.i{1}(find(out{s}.t>=180,1,'first'):find(out{s}.t>=200,1,'first')) *1000);% - mholding_current(t); % + 6 * mean(outleaksub{s}.record{t}.i{1}(find(outleaksub{s}.t>=180,1,'first'):find(outleaksub{s}.t>=200,1,'first')) *1000);
        inewcurr_dend{t,s} =  [out{s}.t';out{s}.record{t}.SEClamp.i{1}' *1000];% - mholding_current(t); % leak subtraction not done by Mongiat + 6 * [0*outleaksub{s}.t';mask.* outleaksub{s}.record{t}.i{1}' *1000];
    end
end
save(fullfile(targetfolder_data,sprintf('Exp_Kir_%s.mat',neuron.experiment)),'mholding_current','neuron','holding_voltage','newcurr_dend','inewcurr_dend','params','vstepsModel','tree','LJP')
