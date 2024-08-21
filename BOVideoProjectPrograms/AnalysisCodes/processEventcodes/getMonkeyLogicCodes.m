%returns also the setting of strobe bit (rising or falling)

function [allcodes_ML,rising,code_start,code_end] = getMonkeyLogicCodes(datapath,filenameML)

try %ML 2
    [b,MLConfig]=mlread([datapath '\' filenameML]);
    code_start=[9];code_end=[18]; %ML2 codes

    if MLConfig.StrobeTrigger==1
        rising=1; %whether strobe bit is rising or falling+
    else
        rising=0;
    end

    %check ML codes
    allcodes_ML=[];
    for i=1:numel(b)
        allcodes_ML=[allcodes_ML;b(i).BehavioralCodes.CodeNumbers b(i).BehavioralCodes.CodeTimes];
    end

catch ME %ML1

    b=bhv_read([datapath '\' filenameML]);
    allcodes_ML=[vertcat(b.CodeNumbers{:}) vertcat(b.CodeTimes{:})];
    code_start=[9 9 9];code_end=[18 18 18]; %ML1 codes
    rising=1; %correct?
end


end