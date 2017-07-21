function handles = aGC_synstim_plot(params,tree,type,ostruct,targetfolder_data,targetfolder_results,showv,c)
if ostruct.newborn
    str = '_newborn';
    str2 = 'bliblablub';
else
   str = '';
   str2 = '_newborn_';
end
if isfield(ostruct,'handles')
    handles = ostruct.handles;
else
    handles = [];
end
if nargin < 8 || isempty(c)
    c = 2;
end
col = colorme(c+1);
folds = dir(targetfolder_data);
folds = {folds.name};
folds = folds(cellfun(@(x) isempty(strfind(x,str2)) & ~isempty(strfind(x,sprintf('Exp_%s_%s%s',type,params.tname,str))),folds));
if isempty(folds)
    errordlg('No experiment found')
    return
else
    answer = listdlg('ListString',folds,'selectionmode','single','promptstring','Please select file to load','listsize',[300,200]);
end
if ~isempty(answer)
    load(fullfile(targetfolder_data,folds{answer}))
else
    return
end
freq = freq; % dämlicher workaround um error zu vermeiden
% freq = [10,20,40,100,150]; % Hz later -> MHz
% dd0 = 0:10:100;  % µm
% dt0 = -50:10:50; % ms
counter = 1;
switch type
    case 'test'
        
    case 'TBS'
        if numel(handles)>=counter && ishandle(handles(counter))
            figure(handles(counter))
        else
            handles(counter) = figure;hold all;
        end
        counter = counter +1;
        for t = 1:numel(tree)-numel(indstim)
            p = plot(out.t,cat(2,out.record{t}.cell.v{recnode{t}}));
            set(p,'Color',tree{t}.col{1})
            p(1).LineStyle = ':';
            p(2).LineStyle = '-.';
            p(3).LineStyle = '-';
            
        end
        legend(p,'Soma','100µm','200µm')
        xlabel('Time [ms]')
        ylabel('Membrane potential [mV]')
        
        if numel(handles)>=counter && ishandle(handles(counter))
            figure(handles(counter))
        else
            handles(counter) = figure;hold all;
        end
        counter = counter +1;
        hold all;
        for t = 1:numel(tree)-numel(indstim)
            p = plot(out.t,out.record{t}.Exp2Syn.i{thesesynids{t}(1)});
            set(p,'Color',tree{t}.col{1})
            
        end
        xlabel('Time [ms]')
        ylabel('Synaptic conductance [nA]')
        
        
        if numel(handles)>=counter && ishandle(handles(counter))
            figure(handles(counter))
        else
            handles(counter) = figure;hold all;
        end
        counter = counter +1;
        hold all;
        for t = 1:numel(tree)-numel(indstim)
            p = plot(out.t,cat(2,out.record{t}.cell.cai{recnode{t}})*1000000);
            set(p,'Color',tree{t}.col{1})
        end
        xlabel('Time [ms]')
        ylabel('Calcium concentration [nM]')
        
    case 'white'
        if numel(handles)>=counter && ishandle(handles(counter))
            figure(handles(counter))
        else
            handles(counter) = figure;hold all;
        end
        counter = counter +1;
        hold all,
        subplot(3,1,1)
        plot(out.t,out.record{t}.cell.v{1})
        subplot(3,1,2)
        f = ones(1,1/neuron{1}.time.dt)*neuron{1}.time.dt/1;
        m_spikeMat = filtfilt(f,1,sum(spikeMat,1));
        plot(tvec,m_spikeMat)
        subplot(3,1,3)
        t2n_plotRaster(spikeMat,tvec)
        
    case 'regular'
        %         handles(end+1) = figure; hold all,
        %         subplot(5,1,1:4)
        %         plot(out.t,out.record{t}.cell.v{1})
        %         subplot(5,1,5)
        % %         line(repmat(1/freq:1/freq:neuron{1}.time.tstop,2,1),repmat([0;1],1,neuron{1}.time.tstop*freq),'color','k','LineWidth',3)
        %         ylim([-0.5 1.5])
        %         set(gca,'ytick',[])
        %         xlim([0,neuron{1}.time.tstop])
        freqmodel = NaN(numel(freq),numel(tree)-1);
        if numel(handles)>=counter && ishandle(handles(counter))
            figure(handles(counter));
        else
            handles(counter) = figure;hold all;
        end
        counter = counter +1;
        fig =[];
        for f = 1:numel(freq)
            subplot(numel(freq),1,f)
            line(repmat(1000/freq(f):1000/freq(f):neuron{1}.time.tstop,2,1),repmat([0;1],1,neuron{1}.time.tstop*freq(f)/1000),'color','k','LineWidth',2)
            for t = 1:numel(tree)-1
                hold all
                %                 plot(out{f}.t,out{f}.record{t}.cell.v{1})
                [~,ind] = findpeaks(out{f}.record{t}.cell.v{1},'MinPeakHeight',0);
                freqmodel(f,t) = numel(ind)/neuron{1}.time.tstop*1000;
                line(repmat(out{f}.t(ind)',2,1),repmat([t;t+1],1,numel(ind)),'color','b','LineWidth',2)
                if freq(f) == 40
                    if numel(handles)>=counter && ishandle(handles(counter))
                        figure(handles(counter));
                    else
                        handles(counter) = figure;hold all;
                    end
                    line(repmat(1000/freq(f):1000/freq(f):neuron{1}.time.tstop,2,1),repmat([0;1],1,neuron{1}.time.tstop*freq(f)/1000),'color','k','LineWidth',2)
                    line(repmat(out{f}.t(ind)',2,1),repmat([t;t+1],1,numel(ind)),'color',col{c},'LineWidth',2)
                    xlim([0 500])
                    figure(handles(counter-1))
                end
            end

            xlim([0,neuron{1}.time.tstop])
            %             ylim([0,numel(tree)+1-1])
        end
        if any(freq == 40)
            figure(handles(counter))
            FontResizer
            FigureResizer(ostruct.figureheight,ostruct.figurewidth)
            counter = counter+1;
        end

        if numel(handles)>=counter && ishandle(handles(counter))
            figure(handles(counter))
        else
            handles(counter) = figure;hold all;
        end
        counter = counter +1;
        errorbar(freq,mean(freqmodel,2),std(freqmodel,[],2))
        line([0, freq(end)],[0, freq(end)],'LineStyle','--','Color',[0.5 0.5 0.5])
        ylim([0, freq(end)])
        xlim([0, freq(end)+5])
        xlabel('Freq_i_n [Hz]')
        ylabel('Freq_o_u_t [Hz]')
    case 'temporal'
        if numel(handles)>=counter && ishandle(handles(counter))
            figure(handles(counter))
        else
            handles(counter) = figure;hold all;
        end
        counter = counter +1;
        freqmodel = NaN(numel(freq),numel(dt0),numel(tree)-2);
        for f = 1:numel(freq)
            for n = 1:numel(dt0)
                ind = (f-1)*numel(dt0) + n;
                subplot(numel(freq),numel(dt0),ind)
                hold all
                line(repmat(1000/freq(f):1000/freq(f):neuron{1}.time.tstop,2,1),repmat([-1;0],1,floor(neuron{1}.time.tstop*freq(f)/1000)),'color','k','LineWidth',2)
                line(repmat(1000/freq(f):1000/freq(f):neuron{1}.time.tstop,2,1)+dt0(n),repmat([-2;-1],1,floor(neuron{1}.time.tstop*freq(f)/1000)),'color','k','LineWidth',2)
                for t = 1:numel(tree)-2
                    if showv
                        plot(out{ind}.t,out{ind}.record{t}.cell.v{1})
                        ylim([-80 numel(tree)])
                    else
                        ylim([0 numel(tree)]-2)
                    end
                    [~,ind2] = findpeaks(out{ind}.record{t}.cell.v{1},'MinPeakHeight',0);
                    freqmodel(f,n,t) = numel(ind2)/neuron{1}.time.tstop*1000;
                    line(repmat(out{ind}.t(ind2)',2,1),repmat([t-1;t],1,numel(ind2)),'color','b','LineWidth',2)
                    if freq(f) == 10
                        if numel(handles)>=counter && ishandle(handles(counter))
                            figure(handles(counter));
                        else
                            handles(counter) = figure;hold all;
                        end
                        if any(dt0(n) == [-15,0])
                            subplot(1,2,find(dt0(n)==[-15,0])), hold all
                            if t == 1
                                line(repmat(1000/freq(f):1000/freq(f):neuron{1}.time.tstop,2,1),repmat([-1;0],1,floor(neuron{1}.time.tstop*freq(f)/1000)),'color','k','LineWidth',2)
                                line(repmat(1000/freq(f):1000/freq(f):neuron{1}.time.tstop,2,1)+dt0(n),repmat([-2;-1],1,floor(neuron{1}.time.tstop*freq(f)/1000)),'color','k','LineWidth',2)
                            end
                            line(repmat(out{ind}.t(ind2)',2,1),repmat([t-1;t],1,numel(ind2)),'color',col{c},'LineWidth',2)
                            xlim([0 500])
                        end
                        figure(handles(counter-1))
                    end
                end
                xlim([0 neuron{1}.time.tstop])
            end
        end
        
        FontResizer
%         FigureResizer(ostruct.figureheight,ostruct.figurewidth,[],ostruct)
        if isfield(ostruct,'savename') && ~isempty(ostruct.savename)
                tprint(fullfile(targetfolder_results,ostruct.savename),'-pdf');
        else
            tprint(fullfile(targetfolder_results,sprintf('TemporalSumRaster_%s_%s%s',type,params.tname,str)),'-pdf');
        end
        if any(freq == 10)
            figure(handles(counter))
            FontResizer
%             FigureResizer(ostruct.figureheight,ostruct.figurewidth)
            counter = counter+1;
        end
        if numel(handles)>=counter && ishandle(handles(counter))
            figure(handles(counter))
        else
            handles(counter) = figure;hold all;
        end
        counter = counter +1;
        col = colorme('dim blue','pink','cyan');
        for f = 1:numel(freq)
            subplot(numel(freq)/2,numel(freq)/2,f)
            hold all,
            title(sprintf('Frequency %g Hz',freq(f)))
%             errorbar(dt0,mean(freqmodel(f,:,:),3)/freq(f),std(freqmodel(f,:,:),[],3)/freq(f),'color',col{f})
            plot(dt0,squeeze(mean(freqmodel(f,:,:),3))/freq(f),'color',col{ostruct.newborn+1})
            ylim([0 1])
            xlim([dt0(1) dt0(numel(dt0))])
            xlabel('\Deltat0 [ms]')
            ylabel('I/O freq. ratio')
        end
        FontResizer
%         FigureResizer(ostruct.figureheight,ostruct.figurewidth,[],ostruct)
        if isfield(ostruct,'savename') && ~isempty(ostruct.savename)
                tprint(fullfile(targetfolder_results,ostruct.savename),'-pdf');
        else
            tprint(fullfile(targetfolder_results,sprintf('TemporalSumPlot_%s_%s%s',type,params.tname,str)),'-pdf');
        end
        if numel(handles)>=counter && ishandle(handles(counter))
            figure(handles(counter))
        else
            handles(counter) = figure;hold all;
        end
        counter = counter +1;
        imagesc(mean(freqmodel,3)./repmat(freq',1,numel(dt0)))
        set(gca,{'CLim','XTick','XTickLabel','YTick','YTickLabel','YDir'},{[0 1],1:numel(dt0),dt0,1:numel(freq),freq,'reverse'})
        xlabel('delta time [ms]')
        ylabel('input frequency [Hz]')
        colorbar
        FontResizer
%         ostruct.image = 1;
%         FigureResizer(5,6,[],ostruct);%ostruct.figureheight,ostruct.figurewidth,[],ostruct)
        if isfield(ostruct,'savename') && ~isempty(ostruct.savename)
                tprint(fullfile(targetfolder_results,ostruct.savename),'-pdf');
        else
            tprint(fullfile(targetfolder_results,sprintf('TemporalSumImag_%s_%s%s',type,params.tname,str)),'-pdf');
        end
    case 'spatial'
        if numel(handles)>=counter && ishandle(handles(counter))
            figure(handles(counter))
        else
            handles(counter) = figure;hold all;
        end
        counter = counter +1;
        freqmodel = NaN(numel(freq),numel(dd0),numel(tree)-1);

        for f = 1:numel(freq)
            for n = 1:numel(dd0)
                ind = (f-1)*numel(dd0) + n;
                subplot(numel(freq),numel(dd0),ind)
                hold all
                line(repmat(1000/freq(f):1000/freq(f):neuron{1}.time.tstop,2,1),repmat([-1;0],1,floor(neuron{1}.time.tstop*freq(f)/1000)),'color','k','LineWidth',2)
                for t = 1:numel(tree)-1
                    if showv
                    plot(out{ind}.t,out{ind}.record{t}.cell.v{1})
                    else
                        ylim([-1 numel(tree)])
                    end
                    [~,ind2] = findpeaks(out{ind}.record{t}.cell.v{1},'MinPeakHeight',0);
                    freqmodel(f,n,t) = numel(ind2)/neuron{1}.time.tstop*1000;
                    line(repmat(out{ind}.t(ind2)',2,1),repmat([t-1;t],1,numel(ind2)),'color','b','LineWidth',2)
                end
                xlim([0 neuron{1}.time.tstop])
            end
        end
        if numel(handles)>=counter && ishandle(handles(counter))
            figure(handles(counter))
        else
            handles(counter) = figure;hold all;
        end
        counter = counter +1;
        hold all
        col = colorme(numel(freq));
        for f = 1:numel(freq)
            subplot(numel(freq),1,f)
%             errorbar(dt0,mean(freqmodel(f,:,:),3)/freq(f),std(freqmodel(f,:,:),[],3)/freq(f),'color',col{f})
            plot(dd0,squeeze(freqmodel(f,:,:))/freq(f),'color',col{f})
            ylim([0 1])
        end
        if numel(handles)>=counter && ishandle(handles(counter))
            figure(handles(counter))
        else
            handles(counter) = figure;hold all;
        end
        counter = counter +1;
        imagesc(mean(freqmodel,3)./repmat(freq',1,numel(dd0)))
        set(gca,{'CLim','XTick','XTickLabel','YTick','YTickLabel'},{[0 1],1:numel(dd0),dd0,1:numel(freq),freq})
        xlabel('delta dist. [µm]')
        ylabel('input frequency [Hz]')
        colorbar
        'g'
    case 'spatial2'
        if numel(handles)>=counter && ishandle(handles(counter))
            figure(handles(counter))
        else
            handles(counter) = figure;hold all;
        end
        counter = counter +1;

        for f = 1:numel(freq)
            for n = 1:numel(dd0)
                ind = (f-1)*numel(dd0) + n;
                ind2 = (f-1)*numel(dd0) + 1;
                subplot(numel(freq),numel(dd0),ind)
                hold all
%                 line(repmat(1000/freq(f):1000/freq(f):neuron{1}.time.tstop,2,1),repmat([-1;0],1,floor(neuron{1}.time.tstop*freq(f)/1000)),'color','k','LineWidth',2)
                for t = 1:numel(tree)-1
                    plot(out{ind}.t,out{ind}.record{t}.cell.v{1})
                end
                xlim([0 neuron{1}.time.tstop])
                ylim([-90 -50])
            end
        end
end
