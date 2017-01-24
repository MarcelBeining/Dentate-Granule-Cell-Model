function [thisdata,steps,rate] = load_ephys(datanum,str,extractkir)

if nargin < 3
    extractkir = 0;
end

folder = 'raw data';  % folder where raw data was placed

switch datanum
    case 1
        load(sprintf('%s%sMongiat_Mature_%s.mat',folder,filesep,str));
        thisdata = data{1};
        if strcmp(str,'CClamp')
            thisdata = thisdata-5;  % don't know why but data is systematically shifted by +5 mV. You see it because it doesnt start at -70 mV
        end
    case 2.21
        load(sprintf('%s%sMongiat_Young_21dpi_%s.mat',folder,filesep,str));
        thisdata = data{1};
    case 2.25
        load(sprintf('%s%sMongiat_Young_25dpi_%s.mat',folder,filesep,str));
        thisdata = data{1};
    case 2.28
        load(sprintf('%s%sMongiat_Young_28dpi_%s.mat',folder,filesep,str));
        thisdata = data{1};
    case 2
        load(sprintf('%s%sMongiat_Young_%s.mat',folder,filesep,str));
        thisdata = data{1};
    case 3
        load(sprintf('%s%sMongiat_BaCl_%s.mat',folder,filesep,str));
        if extractkir && strcmp(str,'VClamp')
            thisdata = data{2}-data{4};  % subtract the Ba insensitive current
        else
            thisdata = data{2};
        end
        
    case 4
        load(sprintf('%s%sMongiat_BaCl_%s.mat',folder,filesep,str));
        if extractkir
            thisdata = data{1}-data{3};  % subtract the Ba insensitive current
        else
            thisdata = data{1};
        end
    case 5
        load(sprintf('%s%sMongiat_BaCl_%s.mat',folder,filesep,str));
        thisdata = data{4};
    case 6
        load(sprintf('%s%sMongiat_BaCl_%s.mat',folder,filesep,str));
        thisdata = data{3};
    case 7
        steps = (50:400)/1000;
        thisdata = NaN(1,1,8);

    case 0
        steps = [];
        thisdata = [];
        rate = [];
    otherwise
        errordlg('datanumber not found')
        return
end

switch str
    case 'VClamp'
        if datanum < 7
            steps = -130:5:-40;
        else
           steps = [];
        end
        
    case 'CClamp'
        if datanum < 7
            steps = (0:5:120)/1000;
        end
end
