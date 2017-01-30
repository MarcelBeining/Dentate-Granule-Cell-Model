function  neuron = changecamp(neuron,cAMP)
if nargin < 2
    cAMP = 1;
end

for t = 1:numel(neuron.mech)
    fields = fieldnames(neuron.mech{t});
    for f1 = 1:numel(fields)
        fields2 = fieldnames(neuron.mech{t}.(fields{f1}));
        if any(strcmp(fields2,'HCN'))
            neuron.mech{t}.(fields{f1}).HCN.cAMP = cAMP;
        end
    end
end


if isfield(neuron,'experiment') && ~isempty(neuron.experiment)
    neuron.experiment = strcat(neuron.experiment,sprintf('_%duM-cAMP',cAMP));
else
    neuron.experiment = sprintf('%duM-cAMP',cAMP);
end

