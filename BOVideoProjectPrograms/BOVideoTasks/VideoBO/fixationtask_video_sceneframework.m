% give names to the TaskObjects defined in the conditions file:
fixation_point = 1;
video = 4;
bgframe = 5;
bg = 6;

mouse_.showcursor(false);

% scene 1: fixation
fix1 = SingleTarget(eye_);
fix1.Target = fixation_point;
fix1.Threshold = 1;
wth1 = WaitThenHold(fix1);
wth1.WaitTime = 5000;
wth1.HoldTime = 300;
scene1 = create_scene(wth1,[fixation_point bgframe bg]);

% scene 2: video
fix2 = SingleTarget(eye_);
fix2.Target = fixation_point;
fix2.Threshold = 1;
wth2 = WaitThenHold(fix2);
wth2.WaitTime = 0;
wth2.HoldTime = 3000;
mg = MovieGraphic(wth2);
scene2 = create_scene(mg,[fixation_point video bgframe bg]);

error_type = 0;
run_scene(scene1,16);
if ~wth1.Success
    if wth1.Waiting
        error_type = 4;  % no fixation
    else
        error_type = 3;  % broke fixation
    end
end

if error_type==0
    run_scene(scene2,23);
    if ~wth2.Success
        error_type = 3;  % broke fixation
    end
end

trialerror(error_type);

if error_type==0
    eventmarker(117); %start reward
    goodmonkey(50,'NumReward', 2, 'PauseTime', 50);
    eventmarker(118); %end reward
end
