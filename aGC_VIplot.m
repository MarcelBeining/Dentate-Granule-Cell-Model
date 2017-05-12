function aGC_VIplot(targetfolder_data,experiment,options)

%  the Kir conductances, we calculate the slope (conductance) at hiperpolarized
% values (-140 to -110mV) and subtracted the slope (conductance) at slightly
% depolarized values (-50 to -70mV). By doing this we calculate only the Kir
% conductance, without the leak conductance at resting potential.
LJP = 11;
if nargin < 3
    options.dataset = 2;
end

if ~isfield(options,'subtract_hv')
    options.subtract_hv = 0;
end
if ~isfield(options,'show')
    options.show = 3;
end


figure;clf;hold all,
colo = colorme({'dim blue','dark red','light blue','red','turquois','orange'});
for o = 1:6
    options.dataset = o;
    switch options.dataset
        case 1
            load('Mongiat_Mature_VClamp.mat');
            exp_vclamp = data{1};
            %             col = colorme('dark blue');
        case 2
            load('Mongiat_Young_VClamp.mat');
            exp_vclamp = data{1};
            %             col = colorme('dark red');
        case 3
            load('Mongiat_BaCl_VClamp.mat');
            exp_vclamp = data{2};
            %             col = colorme('blue');
        case 4
            load('Mongiat_BaCl_VClamp.mat');
            exp_vclamp = data{1};
            %             col = colorme('red');
        case 5
            load('Mongiat_BaCl_VClamp.mat');
            exp_vclamp = data{4};
            %             col = colorme('turquois');
        case 6
            load('Mongiat_BaCl_VClamp.mat');
            exp_vclamp = data{3};
            %             col = colorme('orange');
    end
    tvec = 1/rate:1/rate:size(exp_vclamp,1)/rate;
    
    vstepsreal = vsteps - LJP;
    
    curr_mature = squeeze(mean(exp_vclamp(194*rate+1:204*rate+1,:,:),1));
    basl = squeeze(mean(exp_vclamp(94*rate+1:104*rate+1,:,:),1));
    
    curr_mature = curr_mature-repmat(curr_mature(:,find(vsteps>= -65,1,'first')),1,numel(vsteps));  % -65 war ca das holding pot
    
    spik_ind =cat(2,false(size(curr_mature,1),sum(vsteps<-80)),squeeze(any(exp_vclamp(tvec>0&tvec<204,:,vsteps>=-80) < -300,1))); %index which cells spike
    curr_mature(spik_ind) = NaN;
    
    if any(options.show == [1 3])
        
        mIV = nanmean(curr_mature,1);
        stdIV = nanstd (curr_mature,1);
        hp = patch ([(mIV + stdIV) (fliplr (mIV - stdIV))],[vstepsreal (fliplr (vstepsreal))], colo{o});
        set (hp, 'facealpha', 0.4, 'edgecolor', 'none')
        plot (mIV,vstepsreal, 'Color',colo{o},'LineWidth',3,'LineStyle','.')
    end
    
    xlabel('Measured Current [pA]')
end

for o = 1:6
    switch o
        case 1
            load('Mongiat_Mature_CClamp.mat');
            exp_iclamp = data{1}-5;  % don't know why but data is systematically shifted by +5 mV. You see it because it doesnt start at -70 mV
        case 2
            load('Mongiat_Young_CClamp.mat');
            exp_iclamp = data{1};
        case 3
            load('Mongiat_BaCl_IClamp.mat');
            exp_iclamp = data{2};
        case 4
            load('Mongiat_BaCl_IClamp.mat');
            exp_iclamp = data{1};
        case 5
            load('Mongiat_BaCl_IClamp.mat');
            exp_iclamp = data{4};
        case 6
            load('Mongiat_BaCl_IClamp.mat');
            exp_iclamp = data{3};
    end
    
    % if isfield(ostruct,'variabledt') && ostruct.variabledt == 0
    %     load(expcat(targetfolder_data,'Exp_Spiking',strcat(neuron.experiment,'_fixed-dt')))
    % else
    %     load(expcat(targetfolder_data,'Exp_Spiking',neuron.experiment))
    % end
    tvec = 1/rate:1/rate:size(exp_iclamp,1)/rate;
    %     figure,plot(mean(squeeze(mean(exp_iclamp_mature(tvec<55,:,:),1))))
    vamp = squeeze(mean(exp_iclamp(tvec>205&tvec<255,:,:),1)-0*mean(exp_iclamp(tvec<55,:,:),1));
    vamp(squeeze(any(exp_iclamp(tvec>55&tvec<255,:,:)>0,1))) = NaN;  % delete spikers
    plot(csteps,vamp'-LJP,'color',colo{o})
end

if any(options.show == [2 3])
    if exist(expcat(targetfolder_data,'Exp_Kir',experiment),'file')
        load(expcat(targetfolder_data,'Exp_Kir',experiment))
        
        line(zeros(1,numel(vstepsModel)),vstepsModel,'LineStyle','--','Color',[0.5 0.5 0.5])
        ek = neuron.mech{1}.all.k_ion.ek;
        %         line([-400 100],[ek ek],'LineStyle','--','Color',[0.7 0 0])
        p = plot(steadyStateCurrVec-repmat(steadyStateCurrVec(find(vstepsModel>=-65-LJP,1,'first'),:),size(steadyStateCurrVec,1),1),vstepsModel);
        for t =1:size(steadyStateCurrVec,2)
            set(p(t),'color',tree{t}.col{1})
        end
        ylabel('Holding Voltage [mV] corrected')
    end
    if exist([expcat(targetfolder_data,'Exp_Spiking',experiment),'.mat'],'file')
        
        load([expcat(targetfolder_data,'Exp_Spiking',experiment),'.mat'])
        
        vamp = NaN(numel(cstepsSpikingModel),3);
        for s = 1:numel(cstepsSpikingModel)
            tvec = timeVec{1,s};
            exp_iclamp = cell2mat(voltVec(:,s)');
            if ~any(exp_iclamp(tvec>55&tvec<255,:)>0,1)
                vamp(s,:) = squeeze(mean(exp_iclamp(tvec>205&tvec<255,:),1)-0*mean(exp_iclamp(tvec<55,:),1));
            end
            
        end
        plot(cstepsSpikingModel*1000,vamp,'color',[1 0 1])
    end
end

title('dim blue: mat1, light blue mat2, turqois matBaCl, dark red: young1, light red young2, orange youngBaCl')
FontResizer
% FigureResizer(5,8)
% tprint(fullfile(targetfolder_results,expcat('Fig.2-IV',neuron.experiment)),'-svg');
% tprint(fullfile(targetfolder_results,expcat('Fig.2-IV',neuron.experiment)),'-pdf');
% tprint(fullfile(targetfolder_results,strcat('Fig2-IV',neuron.experiment)),'-png')
% close(fig)
