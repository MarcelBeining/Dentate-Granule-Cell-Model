function aGC_AHP_plot(targetfolder_data,targetfolder_results,neuron)

load(expcat(targetfolder_data,'Exp_msAHP',neuron.experiment))

figure
thisvol = AHPvol_new_curr_dend;
ttthisvol = AHPtvol_new_curr_dend;
% vol_mature = squeeze(mean(mean(exp_vclamp_mature(180*rate+1:200*rate+1,:,:),1),2));
tit = {'Control','XE991 (M Channel blocker)','Apamin (SK blocker)'};
for f=1:size(thisvol,1)
    
    for s = 1:size(thisvol,2)
        subplot(3,1,s)
        hold all
        title(tit{s})
        if ~isempty(ttthisvol{f,s})
            plot(ttthisvol{f,s},squeeze(thisvol{f,s}),'LineWidth',1.5,'Color',tree{f}.col{1})
        end
%         ylim([-70 -30])
%         xlim([0 350])
        line([320,350],[-80 -80],'LineStyle','--','Color','k','LineWidth',5)
        line([600,700],[-80 -80],'Color','k','LineWidth',5)
    end
    
end
linkaxes

for f=1:size(thisvol,1)
    
    for s = 1:size(thisvol,2)
        stimend = find(ttthisvol{f,s} >= 300,1,'first');
        timewindow1 = ttthisvol{f,s} >= 320 & ttthisvol{f,s} <= 350;
        timewindow2 = ttthisvol{f,s} >= 600 & ttthisvol{f,s} <= 700;
        mAHP(s,f) = mean(thisvol{f,s}(timewindow1))-thisvol{f,s}(end);
        sAHP(s,f) = mean(thisvol{f,s}(timewindow2))-thisvol{f,s}(end);
    end
end

display(mAHP)
display(sAHP)