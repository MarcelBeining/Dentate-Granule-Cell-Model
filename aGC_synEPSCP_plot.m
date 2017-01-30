function aGC_synEPSCP_plot(targetfolder_data,targetfolder_results,neuron,ostruct)




if isfield(ostruct,'bablock') && ostruct.bablock
   load(fullfile2(targetfolder_data,'/Exp_synEPSCP_Ba.mat'))
   sim1 = load(fullfile2(targetfolder_data,'/Exp_synEPSCP_CTRL.mat'));
else
    load(fullfile2(targetfolder_data,'/Exp_synEPSCP_CTRL.mat'));
end



figure;hold all;
if isfield(ostruct,'bablock') && ostruct.bablock
    p = plot(mean(abs(sim1.EPSC),2),mean(sim1.EPSP1,2));
    errbar(p,cat(3,std(abs(sim1.EPSP1),[],2),std(sim1.EPSC,[],2)))
end
p = plot(mean(abs(EPSC),2),mean(EPSP1,2));
errbar(p,cat(3,std(abs(EPSP1),[],2),std(EPSC,[],2)))
xlim([0 300]),ylim([0 35])
ylabel('peak EPSP [mV]')
xlabel('peak EPSC [pA]')

tprint(fullfile2(targetfolder_results,expcat('Fig.X-EPSPC',nneuron.experiment)),'-HR-pdf');
