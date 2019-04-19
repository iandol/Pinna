function Pinna_Illusion_grating(radial_speed, rotary_speed)

Screen('Preference', 'SkipSyncTests', 1);
PsychDefaultSetup(2);

screenid = max(Screen('Screens'));
PsychImaging('PrepareConfiguration');
PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');

% [win, winRect] = PsychImaging('OpenWindow', screenid, 0.5, [600, 200, 1300, 900]);
[win, winRect] = PsychImaging('OpenWindow', screenid, 0.5);

ifi = Screen('GetFlipInterval',win);
[win_w, ~] = RectSize(winRect);
[center(1), center(2)] = RectCenter(winRect);

pixelsPerCm = 35;
sizePixel = 10/pixelsPerCm; %mm
viewDistance = 570; %mm
ScreenWidth = round(winRect(3)*sizePixel);
ppd = winRect(3)/2/atand(ScreenWidth/2/viewDistance);
xCen = winRect(3)/2;
yCen = winRect(4)/2;

fix_r = 0.15 ; % radius of fixation point (deg)
fix_cord = [center-fix_r*ppd center+fix_r*ppd];

tilt_angle = 0; %degree
ring_num = 12;
patchNum_perRing = 25;
min_radius = 0.4; %degree
max_radius = (win_w/2+200)/ppd; %degree

offset1 = 0;
ang1 = mod(offset1+linspace(0,360-360/patchNum_perRing,patchNum_perRing),360);
ang1 = deg2rad(ang1);
offset2 = 180/patchNum_perRing;
ang2 = mod(offset2+linspace(0,360-360/patchNum_perRing,patchNum_perRing),360);
ang2 = deg2rad(ang2);

ang_total = repmat([ang1; ang2], ring_num, 1);

sti_duriation = 15; %second

if nargin < 1
    radial_speed = [];
end

if isempty(radial_speed)
    radial_speed = 5; %deg
end

if nargin < 2
    rotary_speed = [];
end

if isempty(rotary_speed)
    rotary_speed = 30; %deg
end

is_mask_outer = 1;
is_mask_inner = 1;
maskinner_radius = 3; %deg
maskouter_radius = 12; %deg

patch_size = [500, 500]; %pixels
sc = 40.0; %Gabor mask parameter
freq = 0.01;
tilt = 0;
contrast = 100.0;
aspectratio = 1.0;
Gabor_parameters = repmat([90, freq, sc, contrast, aspectratio, 0, 0, 0]',1,patchNum_perRing*ring_num);

patch_w = patch_size(1);
patch_h = patch_size(2);
gabortex = CreateProceduralGabor(win, patch_w, patch_h);

try
    white = WhiteIndex(win);
    black = BlackIndex(win);
    gray = (white+black)/2;
    if gray == white
        gray = white/2;
    end
    
    min_radius_pixel = min_radius*ppd;
    sti_radius_pixel = (max_radius - min_radius)*ppd;
    
    texrect = Screen('Rect', gabortex);
    min_patch_size = texrect(3)*0.05;
    max_patch_size = texrect(3)*0.5;
    
    rotary_speed_perFrame = deg2rad(rotary_speed*ifi); %deg/frame
    radial_speed_perFrame = abs(radial_speed)*ppd*ifi; %pixels/frame
    if radial_speed > 0
        onFrames = fix((sti_radius_pixel + sqrt(sti_radius_pixel^2 + 4*sti_radius_pixel*radial_speed_perFrame))/(2*radial_speed_perFrame));
        radial_speed_judgement = 1; %0 means contraction, 1 means expansion, 2 means no radial motion
    elseif radial_speed == 0
        onFrames = fix((sti_radius_pixel + sqrt(sti_radius_pixel^2 + 4*sti_radius_pixel*radial_speed_perFrame))/(2*5*ppd*ifi));
        radial_speed_judgement = 2; %0 means contraction, 1 means expansion, 2 means no radial motion
    elseif radial_speed < 0
        onFrames = fix((sti_radius_pixel + sqrt(sti_radius_pixel^2 + 4*sti_radius_pixel*radial_speed_perFrame))/(2*radial_speed_perFrame));
        radial_speed_judgement = 0; %0 means contraction, 1 means expansion, 2 means no radial motion
    end
    
    maskinner_radius = maskinner_radius * ppd;% convert degrees to pixels
    maskouter_radius = maskouter_radius * ppd;% convert degrees to pixels
    
    if is_mask_outer == 1
        masksize = [winRect(4),winRect(3)];
        mask = ones(masksize(1),masksize(2),2)*gray;
        [x,y] = meshgrid(0:masksize(2)-1,0:masksize(1)-1);
        mask(:, :, 2) = white *(exp(-((x-xCen).^2+(y-yCen).^2)/(maskouter_radius)^2));
        mask(mask>white) = white;
        masktex = Screen('MakeTexture', win, mask);
    end
    if is_mask_inner == 1
        masksize = [winRect(4),winRect(3)];
        mask2 = ones(masksize(1),masksize(2),2)*gray;
        [x,y] = meshgrid(0:masksize(2)-1,0:masksize(1)-1);
        mask2(:, :, 2) = white * (1-exp(-((x-xCen).^2+(y-yCen).^2)/((maskinner_radius)*2)^2));
        mask2(mask2>white) = white;
        masktex2 = Screen('MakeTexture', win, mask2);
    end
    
    n_radial = ones(1,ring_num);
    n_rotary = ones(1,ring_num);
    for i = 1:(ring_num-1)
        n_radial(i+1) = round(i*onFrames/ring_num);
    end
    
    diff_r = 2*sti_radius_pixel/(onFrames)^2;
    diff_size = 2*(max_patch_size-min_patch_size)/(onFrames)^2;
    
    vbl = Screen('Flip',win);
    vblendtime = vbl + sti_duriation;
    
    while vbl < vblendtime
        
        for i = 1:ring_num
            grating_rotAngle = rad2deg(mod(ang_total(i,:)+(n_rotary(i)-1)*rotary_speed_perFrame,2*pi));
            grating_rotAngle_total(:,((i-1)*patchNum_perRing+1):(i*patchNum_perRing)) = grating_rotAngle;
            R = min_radius_pixel + 0.5*diff_r*n_radial(i)^2;
            S = min_patch_size + 0.7*diff_size*n_radial(i)^2;
            
            dstRect(:,((i-1)*patchNum_perRing+1):(i*patchNum_perRing)) = [xCen+R*sind(grating_rotAngle)-S/2; yCen-R*cosd(grating_rotAngle)-S/2;...
                xCen+R*sind(grating_rotAngle)+S/2; yCen-R*cosd(grating_rotAngle)+S/2];
        end
        
        Screen('BlendFunction', win, GL_ONE, GL_ONE);
        Screen('DrawTextures', win, gabortex,[], dstRect, (180-tilt_angle-grating_rotAngle_total), [], [], [], [], kPsychDontDoRotation, Gabor_parameters);
        
        if is_mask_outer == 1
            Screen('BlendFunction', win, GL_ONE_MINUS_SRC_ALPHA,GL_SRC_ALPHA, [1 1 1 0]);
            Screen('DrawTexture', win, masktex);
            Screen('BlendFunction', win, GL_ONE,GL_ZERO, [1 1 1 1]);
        end
        if is_mask_inner == 1
            Screen('BlendFunction', win, GL_ONE_MINUS_SRC_ALPHA,GL_SRC_ALPHA, [1 1 1 0]);
            Screen('DrawTexture', win, masktex2);
            Screen('BlendFunction', win, GL_ONE,GL_ZERO, [1 1 1 1]);
        end
        
        Screen('FillOval', win, [1, 0, 0], fix_cord);
        
        vbl = Screen('Flip', win);
        
        n_rotary = n_rotary + 1;
        if radial_speed_judgement == 1
            n_radial = n_radial + 1;
            n_radial(n_radial == onFrames) = 1;
        elseif radial_speed_judgement == 0
            n_radial = n_radial - 1;
            n_radial(n_radial == 0) = onFrames;
        end
        
        if KbCheck
            break
        end
    end
    sca;
catch
    sca;
end

end