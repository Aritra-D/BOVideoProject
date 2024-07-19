debug = 0;

l=load('stimparams_BO.mat','stimdur','scenespertrial','blankdur');

% give names to the TaskObjects defined in the conditions file:
% fiximage_fin = 1;
fiximage = 1;
stimfigs = 2:2+l.scenespertrial-1;
blankscreen = 2+l.scenespertrial;
bgscreen = 2+l.scenespertrial+1;

% define time intervals (in ms):
prestimtime = 300;
targetappeartime = l.stimdur; %O Herron: analysis window 100-500 ms
blanktime = l.blankdur;

maxwaittime = 10000/1e3;
waittimeforsacc = 3000;

if debug
    waittimeforsacc = 500;
end

%reward parameters
rewardtime = 50;
rewardnum = 2;
rewardpause = 100;

% fixation windows:
fix_radius = 1; %1 for V4 

%acquire fix
toggleobject([fiximage blankscreen bgscreen],'eventmarker',23,'status','on');

%allow monkey to self-initiate trial
fixationpointheld=0;
t=tic;
while ~fixationpointheld
    [ontarget,rt]= eyejoytrack('acquirefix', fiximage, fix_radius, 1000);
    if(ontarget||debug)
        eventmarker(1); %acquired fixation
        [ontarget]= eyejoytrack('holdfix', fiximage, fix_radius, prestimtime);
        if(ontarget||debug)
            if debug
               idle(prestimtime); 
            end
            eventmarker(6);
            fixationpointheld=1;
        end
    end
    if(toc(t)>maxwaittime && fixationpointheld~=1)
        disp('too long')
        eventmarker(13);
        trialerror(2);
      
        toggleobject(fiximage,'eventmarker',26,'status','off');
        return;
    end
end

for i=1:numel(stimfigs) %play all square scenes in the trial with blank interval
    %scene on
    toggleobject([stimfigs(i)],'eventmarker',(TrialRecord.CurrentBlock)+100,'status','on');
    [ontarget]= eyejoytrack('holdfix', fiximage, fix_radius, targetappeartime);
    if(debug)
        idle(targetappeartime);
        ontarget=1;
    end
    if ~ontarget
        disp('here')
        eventmarker(14);
        trialerror(3);

        toggleobject([stimfigs fiximage],'eventmarker',26,'status','off');
        return;
    end

    %scene off (blank screen was already on)
    toggleobject([stimfigs(i)],'eventmarker',26);
    [ontarget]= eyejoytrack('holdfix', fiximage, fix_radius, blanktime);
    if(debug)
        idle(blanktime);
        ontarget=1;
    end
    if ~ontarget
        disp('here')
        eventmarker(14);
        trialerror(3);

        toggleobject([stimfigs fiximage],'eventmarker',26,'status','off');
        return;
    end
end

if ~ontarget
   eventmarker(14);
   trialerror(3);
   toggleobject([fiximage stimfig  stimambig],'eventmarker',26,'status','off');
   return;
else 
   toggleobject([fiximage],'eventmarker',6,'status','off');
   trialerror(0);
   goodmonkey(rewardtime,'NumReward',rewardnum, 'PauseTime', rewardpause);
   eventmarker(117);
end