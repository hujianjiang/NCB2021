clear
clc
close all
tic


filepath = 'Y:\Optogenetics\HT1080PARac1T17N\analysis\windowing\';
InitialArray = 1:33;
EXPname = 'HT1080PARac1T17Nin1%FBSon2.3kPAA';
savefolder = 'results20210203_FORCE2\';
mkdir([filepath(InitialArray) savefolder]);

TimePoints = 59 : -1 : 56;
TimePointsSelected = [59,58];
PlotRange = 50 : 87;
SPEEDlimit = 1;
SmoothFactor = 0.01;

timeinterval = 5;
distance = 12;
FORCEdepth = [2,2];
valuerange = 5;
protArray = [];
retrArray = [];

protSamples{5} = load([filepath EXPname sprintf('%03.f',5) '\WindowingPackage\protrusion_samples\protrusion_samples.mat']);
[m,n] = size(protSamples{5}.protSamples.avgNormal);
FORCEprotrusion(1,:) = nan(1,2*n);
FORCEretraction(1,:) = nan(1,2*n);
SPEEDprotrusion(1,:) = nan(1,2*n);
SPEEDretraction(1,:) = nan(1,2*n);
Accelprotrusion(1,:) = nan(1,2*n);
Accelretraction(1,:) = nan(1,2*n);
InteDprotrusion(1,:) = nan(1,2*n);
InteDretraction(1,:) = nan(1,2*n);
lengthafterretraction(1) = NaN;
lengthafterprotrusion(1) = NaN;
FORCEprotrusionMaxV(1,:) = nan(1,2*n);
FORCEretractionMaxV(1,:) = nan(1,2*n);
SPEEDprotrusionMaxV(1,:) = nan(1,2*n);
SPEEDretractionMaxV(1,:) = nan(1,2*n);
AccelprotrusionMaxV(1,:) = nan(1,2*n);
AccelretractionMaxV(1,:) = nan(1,2*n);
InteDprotrusionMaxV(1,:) = nan(1,2*n);
InteDretractionMaxV(1,:) = nan(1,2*n);

for i = [6, 8:38, 40:49, 52:78] 

    % load original data
    protSamples{i} = load([filepath EXPname sprintf('%03.f',i) '\WindowingPackage\protrusion_samples\protrusion_samples.mat']);
    TfSamples{i} = load([filepath EXPname sprintf('%03.f',i) '\WindowingPackage\window_sampling\Raw images - channel 1.mat']);
    for j1 = 61 : 84
        windows{i,j1} = load([filepath EXPname sprintf('%03.f',i) '\WindowingPackage\windows\windows_frame__frame_' sprintf('%03.f',j1) '.mat']);
    end
    ROI{i} = importdata([filepath(InitialArray) 'ROI_aligned\' EXPname sprintf('%03.f',i) '_ROIaligned.mat']);
    
    % Align the mismached timeframes by window sampling
    FORCEchosen = squeeze(TfSamples{i}.samples.avg(:,1,:));
    [m3,n3] = size(FORCEchosen);
    for i4 = 2 : n3
        FORCEchosena = FORCEchosen(:,i4-1);
        FORCEchosenb = FORCEchosen(:,i4);
        FORCEchosenaTEMP = repnan(FORCEchosena,'pchip');
        FORCEchosenbTEMP = repnan(FORCEchosenb,'pchip');
        [r,lag] = xcorr(FORCEchosenaTEMP,FORCEchosenbTEMP);
        
        shift{i4} = lag(r==max(r));
        if find(length(shift{i4}))
            FORCEchosen(:,i4) = circshift(FORCEchosenb,shift{i4}(1));
        else
            FORCEchosena = FORCEchosen(:,i4-2);
            FORCEchosenb = FORCEchosen(:,i4);
            FORCEchosenaTEMP = repnan(FORCEchosena,'pchip');
            FORCEchosenbTEMP = repnan(FORCEchosenb,'pchip');
            [r,lag] = xcorr(FORCEchosenaTEMP,FORCEchosenbTEMP);
            
            shift{i4} = lag(r==max(r));
            FORCEchosen(:,i4) = circshift(FORCEchosenb,shift{i4});
        end
        
        protSamples{i}.protSamples.avgNormal(:,i4) = circshift(protSamples{i}.protSamples.avgNormal(:,i4),shift{i4});
        for i6 = 1 : size(TfSamples{i}.samples.avg,2)
            TfSamples{i}.samples.avg(:,i6,i4) = circshift(TfSamples{i}.samples.avg(:,i6,i4),shift{i4});
        end
    end
    
    
    % generate the integrated protrusion map
    [m2,n2] = size(protSamples{i}.protSamples.avgNormal);
    for i0 = 1 : n2-1
        protIntegrated{i}(:,i0) = nansum(protSamples{i}.protSamples.avgNormal(:,1:i0),2);
    end
    
    
    % find the windows within the PA-ROI region
    WindowsROI{i}=[];
    for j2 = 1 : length(windows{i, 61}.windows)
        if ~isempty(windows{i, 61}.windows{1,j2})
            clearvars temp1
            temp1 = [windows{i, 61}.windows{1, j2}{1, 1}{1, 1}  windows{i, 61}.windows{1, j2}{1, 1}{1, 2} windows{i, 61}.windows{1, j2}{1, 1}{1, 3} windows{i, 61}.windows{1, j2}{1, 1}{1, 4}];
            if find(~ROI{1,i}{61}(round(temp1(2,:)),round(temp1(1,:))))
                continue
            else
                WindowsROI{i} = [WindowsROI{i} j2];
            end
        end
    end
    
    %     % remove the windows move out of the ROI region during PA
        if ~isempty(WindowsROI{i})
            tempPos = [];
            for j2 = 1 : length(WindowsROI{i})
                for j21 = 61 : 84
                    if WindowsROI{i}(j2) > length(windows{i, j21}.windows) || isempty(windows{i, j21}.windows{1, WindowsROI{i}(j2)})
                        tempPos = [tempPos, j2];
                        break
                    end
                    clearvars temp1
                    temp1 = [windows{i, j21}.windows{1, WindowsROI{i}(j2)}{1, 1}{1, 1}  windows{i, j21}.windows{1, WindowsROI{i}(j2)}{1, 1}{1, 2} windows{i, j21}.windows{1, WindowsROI{i}(j2)}{1, 1}{1, 3} windows{i, j21}.windows{1, WindowsROI{i}(j2)}{1, 1}{1, 4}];
                    if find(ROI{1,i}{j21}(round(temp1(2,:)),round(temp1(1,:))))
                        continue
                    else
                        tempPos = [tempPos, j2];
                        break
                    end
                end
            end
            if ~isempty(tempPos)
                WindowsROI{i}(tempPos) = [];
            end
        end
    
    
    
    % check the time range of the PA-ROI windows before PA
    for j3 = WindowsROI{i}
        temp = fillmissing(protIntegrated{i}(j3, :), 'movmedian',20);
        [pp,p] = csaps(1:n-1 , temp, SmoothFactor);
        % get smoothed membrane dynamic & speed
        protIntegratedSmoothed{i}(j3,:) = pp.coefs(:,4);
        protSamplesSmoothed{i}(j3,1) = pp.coefs(1,4);
        for k1 = 2 : n-2
            protSamplesSmoothed{i}(j3,k1) = pp.coefs(k1,4) - pp.coefs(k1-1,4);
        end
        % get smoothed & unsmoothed membrane speed acceleration
        for k1 = 1 : n-1
            protAcceleration{i}(j3,k1) = protSamples{i}.protSamples.avgNormal(j3,k1+1) - protSamples{i}.protSamples.avgNormal(j3,k1);
        end
            temp2 = fillmissing(protAcceleration{i}(j3,:),'movmedian',20);
            [pp1,p1] = csaps(1:n-1,temp2,SmoothFactor);
            if length(pp1.coefs(:,4))==n-2
                protAccelerationSmoothed{i}(j3,1:n-2) = pp1.coefs(:,4);
            else
                protAccelerationSmoothed{i}(j3,1:n-2) = NaN;
            end
        
        
        
        % get smoothed traction force
        temp = nanmean([squeeze(TfSamples{i}.samples.avg(j3,FORCEdepth(1),:)),squeeze(TfSamples{i}.samples.avg(j3,FORCEdepth(2),:))],2);
        tempTF = fillmissing(temp,'movmedian',20);
        
        TfSamplesSmoothed{i}(j3,1:n-1) = tempTF(1:n-1);

        
        if size(pp.coefs,2) == 4
            [peakval, peakloc, peakwidth, peakprominence] = findpeaks(pp.coefs(1 : 61,4));
            
            peaklocWithEnds{i}{j3} = [1; peakloc; 61];
            peakvalWithEnds{i}{j3} = [pp.coefs(1,4); peakval; pp.coefs(61,4);];
            
            % find local valley between each peak pairs
            for k1 = 1 : length(peaklocWithEnds{i}{j3})-1
                [minval, minidx] = min(pp.coefs(peaklocWithEnds{i}{j3}(k1):peaklocWithEnds{i}{j3}(k1+1),4));
                
                valleyloc{i}{j3}(k1) = peaklocWithEnds{i}{j3}(k1) + minidx - 1;
                valleyval{i}{j3}(k1) = minval;
                if valleyloc{i}{j3}(k1) == 61
                    valleyloc{i}{j3}(k1) = [];
                    valleyval{i}{j3}(k1) = [];
                end
                
            end
            
            % remove the short term switches
            peakvalfilted{i}{j3} = peakval;
            peaklocfilted{i}{j3} = peakloc;
            for k2 = 1 : n2-distance
                temploc = find(peakloc>k2 & peakloc<=k2+distance);
                temptemploc = find(peakval(temploc) ~= max(peakval(temploc)));
                peaklocfilted{i}{j3}(temploc(temptemploc)) = NaN;
                peakvalfilted{i}{j3}(temploc(temptemploc)) = NaN;
            end
            
            
            
            % get the protrusion & retraction arrays
            peakposition = find(~isnan(peaklocfilted{i}{j3}));
            if length(peakposition) > 1
                for i3 = 1 : length(peakposition)-1
                    tempvalleyloc = find(valleyloc{i}{j3}>peaklocfilted{i}{j3}(peakposition(i3)) & valleyloc{i}{j3}<peaklocfilted{i}{j3}(peakposition(i3+1)));
                    % add retraction array if its time length is larger than given distance and its range is larger than given scale.
                    if peakposition(i3) > 1 
                        
                        FORCEretraction = [FORCEretraction ; nan(1,2*n)];
                        InteDretraction = [InteDretraction; nan(1,2*n)];
                        SPEEDretraction = [SPEEDretraction; nan(1,2*n)];
                        Accelretraction = [Accelretraction; nan(1,2*n)];
                        lengthafterretraction = [lengthafterretraction; NaN];
                        
                        lengthafterretraction(end) = valleyloc{i}{j3}(tempvalleyloc(1))-peaklocfilted{i}{j3}(peakposition(i3))+1;
                        
                        FORCEretraction(end,n+1:n+lengthafterretraction(end)) = TfSamplesSmoothed{i}(j3,peaklocfilted{i}{j3}(peakposition(i3)):valleyloc{i}{j3}(tempvalleyloc(1)));
                        InteDretraction(end,n+1:n+lengthafterretraction(end)) = pp.coefs(peaklocfilted{i}{j3}(peakposition(i3)):valleyloc{i}{j3}(tempvalleyloc(1)),4);
                        SPEEDretraction(end,n+1:n+lengthafterretraction(end)) = protSamplesSmoothed{i}(j3,peaklocfilted{i}{j3}(peakposition(i3)):valleyloc{i}{j3}(tempvalleyloc(1)));
                        Accelretraction(end,n+1:n+lengthafterretraction(end)) = protAccelerationSmoothed{i}(j3,peaklocfilted{i}{j3}(peakposition(i3)):valleyloc{i}{j3}(tempvalleyloc(1)));
                        
                        clearvars retractionV pp1 p1 m1 n1
                        retractionV = protSamplesSmoothed{i}(j3,peaklocfilted{i}{j3}(peakposition(i3)):valleyloc{i}{j3}(tempvalleyloc(1)));
                        [pp1,p1] = csaps(1:length(retractionV),retractionV,SmoothFactor);
                        if size(pp1.coefs,2) == 4
                            [m1,n1] = find(pp1.coefs(:,4) == min(pp1.coefs(:,4)));
                            
                            
                            FORCEretractionMaxV = [FORCEretractionMaxV ; nan(1,2*n)];
                            InteDretractionMaxV = [InteDretractionMaxV ; nan(1,2*n)];
                            SPEEDretractionMaxV = [SPEEDretractionMaxV ; nan(1,2*n)];
                            AccelretractionMaxV = [AccelretractionMaxV ; nan(1,2*n)];
                            
                            
                            FORCEretractionMaxV(end,n-m1-peaklocfilted{i}{j3}(peakposition(i3))+2:n) = TfSamplesSmoothed{i}(j3,1:peaklocfilted{i}{j3}(peakposition(i3))+m1-1);
                            FORCEretractionMaxV(end,n+1:n+62-m1-peaklocfilted{i}{j3}(peakposition(i3))) = TfSamplesSmoothed{i}(j3,peaklocfilted{i}{j3}(peakposition(i3))+m1:61);
                            SPEEDretractionMaxV(end,n-m1-peaklocfilted{i}{j3}(peakposition(i3))+2:n) = protSamplesSmoothed{i}(j3,1:peaklocfilted{i}{j3}(peakposition(i3))+m1-1);
                            SPEEDretractionMaxV(end,n+1:n+62-m1-peaklocfilted{i}{j3}(peakposition(i3))) = protSamplesSmoothed{i}(j3,peaklocfilted{i}{j3}(peakposition(i3))+m1:61);
                            AccelretractionMaxV(end,n-m1-peaklocfilted{i}{j3}(peakposition(i3))+2:n) = protAccelerationSmoothed{i}(j3,1:peaklocfilted{i}{j3}(peakposition(i3))+m1-1);
                            AccelretractionMaxV(end,n+1:n+62-m1-peaklocfilted{i}{j3}(peakposition(i3))) = protAccelerationSmoothed{i}(j3,peaklocfilted{i}{j3}(peakposition(i3))+m1:61);
                            
                            InteDretractionMaxV(end,n-m1-peaklocfilted{i}{j3}(peakposition(i3))+2:n) = pp.coefs(1:peaklocfilted{i}{j3}(peakposition(i3))+m1-1,4);
                            InteDretractionMaxV(end,n+1:n+62-m1-peaklocfilted{i}{j3}(peakposition(i3))) = pp.coefs(peaklocfilted{i}{j3}(peakposition(i3))+m1:61,4);
                        end
                        
                        %                                                 end
                    end
                    
                    % add protrusion array if its length is larger than given distance.
                    if peakposition(i3+1) <=61 
                        FORCEprotrusion = [FORCEprotrusion ; nan(1,2*n)];
                        InteDprotrusion = [InteDprotrusion; nan(1,2*n)];
                        SPEEDprotrusion = [SPEEDprotrusion; nan(1,2*n)];
                        Accelprotrusion = [Accelprotrusion; nan(1,2*n)];
                        lengthafterprotrusion = [lengthafterprotrusion; NaN];
                        
                        lengthafterprotrusion(end) = peaklocfilted{i}{j3}(peakposition(i3+1))-valleyloc{i}{j3}(tempvalleyloc(1))+1;
                        
                        
                        FORCEprotrusion(end,n+1:n+lengthafterprotrusion(end)) = TfSamplesSmoothed{i}(j3,valleyloc{i}{j3}(tempvalleyloc(1)):peaklocfilted{i}{j3}(peakposition(i3+1)));
                        InteDprotrusion(end,n+1:n+lengthafterprotrusion(end)) = pp.coefs(valleyloc{i}{j3}(tempvalleyloc(1)):peaklocfilted{i}{j3}(peakposition(i3+1)),4);
                        SPEEDprotrusion(end,n+1:n+lengthafterprotrusion(end)) = protSamplesSmoothed{i}(j3,valleyloc{i}{j3}(tempvalleyloc(1)):peaklocfilted{i}{j3}(peakposition(i3+1)));
                        Accelprotrusion(end,n+1:n+lengthafterprotrusion(end)) = protAccelerationSmoothed{i}(j3,valleyloc{i}{j3}(tempvalleyloc(1)):peaklocfilted{i}{j3}(peakposition(i3+1)));
                        
                        
                        clearvars protrusionV pp1 p1 m1 n1
                        %                         protrusionV = protSamples{i}.protSamples.avgNormal(j3,valleyloc{i}{j3}(tempvalleyloc(1)):peaklocfilted{i}{j3}(peakposition(i3+1)));
                        protrusionV = protSamplesSmoothed{i}(j3,valleyloc{i}{j3}(tempvalleyloc(1)):peaklocfilted{i}{j3}(peakposition(i3+1)));
                        if length(find(~isnan(protrusionV))) > 1
                            [pp1,p1] = csaps(1:length(protrusionV),protrusionV,SmoothFactor);
                            if size(pp1.coefs,2) == 4
                                FORCEprotrusionMaxV = [FORCEprotrusionMaxV ; nan(1,2*n)];
                                InteDprotrusionMaxV = [InteDprotrusionMaxV ; nan(1,2*n)];
                                SPEEDprotrusionMaxV = [SPEEDprotrusionMaxV ; nan(1,2*n)];
                                AccelprotrusionMaxV = [AccelprotrusionMaxV ; nan(1,2*n)];
                                
                                [m1,n1] = find(pp1.coefs(:,4) == max(pp1.coefs(:,4)));
                                
                                FORCEprotrusionMaxV(end,n-valleyloc{i}{j3}(tempvalleyloc(1))-m1+2:n) = TfSamplesSmoothed{i}(j3,1:valleyloc{i}{j3}(tempvalleyloc(1))+m1-1);
                                FORCEprotrusionMaxV(end,n+1:n+62-valleyloc{i}{j3}(tempvalleyloc(1))-m1) = TfSamplesSmoothed{i}(j3,valleyloc{i}{j3}(tempvalleyloc(1))+m1:61);
                                SPEEDprotrusionMaxV(end,n-valleyloc{i}{j3}(tempvalleyloc(1))-m1+2:n) = protSamplesSmoothed{i}(j3,1:valleyloc{i}{j3}(tempvalleyloc(1))+m1-1);
                                SPEEDprotrusionMaxV(end,n+1:n+62-valleyloc{i}{j3}(tempvalleyloc(1))-m1) =protSamplesSmoothed{i}(j3,valleyloc{i}{j3}(tempvalleyloc(1))+m1:61);
                                AccelprotrusionMaxV(end,n-valleyloc{i}{j3}(tempvalleyloc(1))-m1+2:n) = protAccelerationSmoothed{i}(j3,1:valleyloc{i}{j3}(tempvalleyloc(1))+m1-1);
                                AccelprotrusionMaxV(end,n+1:n+62-valleyloc{i}{j3}(tempvalleyloc(1))-m1) =protAccelerationSmoothed{i}(j3,valleyloc{i}{j3}(tempvalleyloc(1))+m1:61);
                                
                                InteDprotrusionMaxV(end,n-valleyloc{i}{j3}(tempvalleyloc(1))-m1+2:n) = pp.coefs(1:valleyloc{i}{j3}(tempvalleyloc(1))+m1-1,4);
                                InteDprotrusionMaxV(end,n+1:n+62-valleyloc{i}{j3}(tempvalleyloc(1))-m1) = pp.coefs(valleyloc{i}{j3}(tempvalleyloc(1))+m1:61,4);
                                
                                
                            end
                        end
                        
                        
                    end
                    if valleyloc{i}{j3}(end) == 61
                        valleyloc{i}{j3}(end) = [];
                        valleyval{i}{j3}(end) = [];
                    end
                    
                end
            end
            
        end
        protTemp = find(valleyloc{i}{j3}<=61);
        retrTemp = find(peaklocfilted{i}{j3}<=61);
        if isempty(retrTemp) && isempty(protTemp)
            continue
        else if isempty(retrTemp)
                protArray = [protArray; i j3 valleyloc{i}{j3}(protTemp(end))];
            else if isempty(protTemp)
                    retrArray = [retrArray; i j3 peaklocfilted{i}{j3}(retrTemp(end))];
                else if valleyloc{i}{j3}(protTemp(end)) > peaklocfilted{i}{j3}(retrTemp(end))
                        protArray = [protArray; i j3 valleyloc{i}{j3}(protTemp(end))];
                    else if valleyloc{i}{j3}(protTemp(end)) < peaklocfilted{i}{j3}(retrTemp(end))
                            retrArray = [retrArray; i j3 peaklocfilted{i}{j3}(retrTemp(end))];
                        end
                    end
                end
            end
        end
        
    end
    
        temp1 = find(protArray(:,1)==i);
        temp2 = find(retrArray(:,1)==i);
    
    
    [i toc/60]
end





% data select & normalize
for i1 = TimePoints
    protSegments{i1} = [];
    protSegmentsSpeed{i1} = [];
    protSegmentsTf{i1} = [];
    protSegmentsAccel{i1} = [];
    
    temp2 = find(protArray(:,3) <= i1 & protArray(:,3) >= i1);
    if ~isempty(temp2)
        for j4 = temp2'
            protSegments{i1} = [protSegments{i1}; protIntegratedSmoothed{1, protArray(j4,1)}(protArray(j4,2),:)];
            protSegmentsSpeed{i1} = [protSegmentsSpeed{i1}; protSamplesSmoothed{1, protArray(j4,1)}(protArray(j4,2),:)];
            tempTf = TfSamplesSmoothed{1, protArray(j4,1)};
            protSegmentsTf{i1} = [protSegmentsTf{i1}; tempTf(protArray(j4,2),:)];
            protSegmentsAccel{i1} = [protSegmentsAccel{i1}; protAccelerationSmoothed{1, protArray(j4,1)}(protArray(j4,2),:)];

            
        end
    end
    
    
    retrSegments{i1} = [];
    retrSegmentsSpeed{i1} = [];
    retrSegmentsTf{i1} = [];
    retrSegmentsAccel{i1} = [];
    
    temp3 = find(retrArray(:,3) <= i1 & retrArray(:,3) >= i1);
    
    if ~isempty(temp3)
        for j5 = temp3'
            retrSegments{i1} = [retrSegments{i1}; protIntegratedSmoothed{1, retrArray(j5,1)}(retrArray(j5,2),:)];
            retrSegmentsSpeed{i1} = [retrSegmentsSpeed{i1}; protSamplesSmoothed{1, retrArray(j5,1)}(retrArray(j5,2),:)];
            tempTf = TfSamplesSmoothed{1, retrArray(j5,1)};
            retrSegmentsTf{i1} = [retrSegmentsTf{i1}; tempTf(retrArray(j5,2),:)];
            retrSegmentsAccel{i1} = [retrSegmentsAccel{i1}; protAccelerationSmoothed{1, retrArray(j5,1)}(retrArray(j5,2),:)];
            
        end
    end
    
end


% Data select
InteDprotrusionSelect = InteDprotrusion(:,n+1:n+length(PlotRange));
SPEEDprotrusionSelect = SPEEDprotrusion(:,n+1:n+length(PlotRange));
AccelprotrusionSelect = Accelprotrusion(:,n+1:n+length(PlotRange));
FORCEprotrusionSelect = FORCEprotrusion(:,n+1:n+length(PlotRange));
temp5 = [];
for i5 = 1 : size(SPEEDprotrusionSelect,1)
    if ~isempty(find((abs(SPEEDprotrusionSelect(i5,:))>SPEEDlimit)))
        temp5 = [temp5,i5];
    end
end
if ~isempty(temp5)
    InteDprotrusionSelect(temp5,:) = [];
    SPEEDprotrusionSelect(temp5,:) = [];
    AccelprotrusionSelect(temp5,:) = [];
    FORCEprotrusionSelect(temp5,:) = [];
end



InteDretractionSelect = -InteDretraction(:,n+1:n+length(PlotRange));
SPEEDretractionSelect = -SPEEDretraction(:,n+1:n+length(PlotRange));
AccelretractionSelect = Accelretraction(:,n+1:n+length(PlotRange));
FORCEretractionSelect = FORCEretraction(:,n+1:n+length(PlotRange));
temp5 = [];
for i5 = 1 : size(SPEEDretractionSelect,1)
    if ~isempty(find((abs(SPEEDretractionSelect(i5,:))>SPEEDlimit)))
        temp5 = [temp5,i5];
    end
end
if ~isempty(temp5)
    InteDretractionSelect(temp5,:) = [];
    SPEEDretractionSelect(temp5,:) = [];
    AccelretractionSelect(temp5,:) = [];
    FORCEretractionSelect(temp5,:) = [];
end

protN1 = find(nanmean(InteDprotrusionMaxV(:,n-10:n))==min(nanmean(InteDprotrusionMaxV(:,n-10:n))));
InteDprotrusionMaxVSelect = InteDprotrusionMaxV(:,n-10+protN1 : n+protN1+21);
SPEEDprotrusionMaxVSelect = SPEEDprotrusionMaxV(:,n-10+protN1 : n+protN1+21);
AccelprotrusionMaxVSelect = AccelprotrusionMaxV(:,n-10+protN1 : n+protN1+21);
FORCEprotrusionMaxVSelect = FORCEprotrusionMaxV(:,n-10+protN1 : n+protN1+21);
temp5 = [];
for i5 = 1 : size(SPEEDprotrusionMaxVSelect,1)
    if ~isempty(find((abs(SPEEDprotrusionMaxVSelect(i5,:))>SPEEDlimit)))
        temp5 = [temp5,i5];
    end
end
if ~isempty(temp5)
    InteDprotrusionMaxVSelect(temp5,:) = [];
    SPEEDprotrusionMaxVSelect(temp5,:) = [];
    AccelprotrusionMaxVSelect(temp5,:) = [];
    FORCEprotrusionMaxVSelect(temp5,:) = [];
end

retrN1 = find(nanmean(InteDretractionMaxV(:,n-10:n))==max(nanmean(InteDretractionMaxV(:,n-10:n))));
InteDretractionMaxVSelect = -InteDretractionMaxV(:,n-11+retrN1 : n+retrN1+20);
SPEEDretractionMaxVSelect = -SPEEDretractionMaxV(:,n-11+retrN1 : n+retrN1+20);
AccelretractionMaxVSelect = AccelretractionMaxV(:,n-11+retrN1 : n+retrN1+20);
FORCEretractionMaxVSelect = FORCEretractionMaxV(:,n-11+retrN1 : n+retrN1+20);
temp5 = [];
for i5 = 1 : size(SPEEDretractionMaxVSelect,1)
    if ~isempty(find((abs(SPEEDretractionMaxVSelect(i5,:))>SPEEDlimit)))
        temp5 = [temp5,i5];
    end
end
if ~isempty(temp5)
    InteDretractionMaxVSelect(temp5,:) = [];
    SPEEDretractionMaxVSelect(temp5,:) = [];
    AccelretractionMaxVSelect(temp5,:) = [];
    FORCEretractionMaxVSelect(temp5,:) = [];
end



for i1 =  TimePoints
    protSegmentsSelect{i1} = protSegments{i1}(:,i1 : i1+length(PlotRange));
    protSegmentsSpeedSelect{i1} = protSegmentsSpeed{i1}(:,i1 : i1+length(PlotRange));
    protSegmentsAccelSelect{i1} = protSegmentsAccel{i1}(:,i1 : i1+length(PlotRange));
    protSegmentsTfSelect{i1} = protSegmentsTf{i1}(:,i1 : i1+length(PlotRange));
    
    temp5 = [];
    for i5 = 1 : size(protSegmentsSpeedSelect{i1},1)
        if ~isempty(find((abs(protSegmentsSpeedSelect{i1}(i5,:))>SPEEDlimit)))
            temp5 = [temp5,i5];
        end
    end
    if ~isempty(temp5)
        protSegmentsSelect{i1}(temp5,:) = [];
        protSegmentsSpeedSelect{i1}(temp5,:) = [];
        protSegmentsAccelSelect{i1}(temp5,:) = [];
        protSegmentsTfSelect{i1}(temp5,:) = [];
    end
    
    
    retrSegmentsSelect{i1} = -retrSegments{i1}(:,i1 : i1+length(PlotRange));
    retrSegmentsSpeedSelect{i1} = -retrSegmentsSpeed{i1}(:,i1 : i1+length(PlotRange));
    retrSegmentsAccelSelect{i1} = retrSegmentsAccel{i1}(:,i1 : i1+length(PlotRange));
    retrSegmentsTfSelect{i1} = retrSegmentsTf{i1}(:,i1 : i1+length(PlotRange));
    temp5 = [];
    for i5 = 1 : size(retrSegmentsSpeedSelect{i1},1)
        if ~isempty(find((abs(retrSegmentsSpeedSelect{i1}(i5,:))>SPEEDlimit)))
            temp5 = [temp5,i5];
        end
    end
    if ~isempty(temp5)
        retrSegmentsSelect{i1}(temp5,:) = [];
        retrSegmentsSpeedSelect{i1}(temp5,:) = [];
        retrSegmentsAccelSelect{i1}(temp5,:) = [];
        retrSegmentsTfSelect{i1}(temp5,:) = [];
    end
    
end



% Data normalize
InteDprotrusionNorm = normalize(InteDprotrusionSelect,2);
SPEEDprotrusionNorm = normalize(SPEEDprotrusionSelect,2);
AccelprotrusionNorm = normalize(AccelprotrusionSelect,2);
FORCEprotrusionNorm = normalize(FORCEprotrusionSelect,2);

InteDretractionNorm = normalize(InteDretractionSelect,2);
SPEEDretractionNorm = normalize(SPEEDretractionSelect,2);
AccelretractionNorm = normalize(AccelretractionSelect,2);
FORCEretractionNorm = normalize(FORCEretractionSelect,2);

InteDprotrusionMaxVNorm = normalize(InteDprotrusionMaxVSelect,2);
SPEEDprotrusionMaxVNorm = normalize(SPEEDprotrusionMaxVSelect,2);
AccelprotrusionMaxVNorm = normalize(AccelprotrusionMaxVSelect,2);
FORCEprotrusionMaxVNorm = normalize(FORCEprotrusionMaxVSelect,2);

InteDretractionMaxVNorm = normalize(InteDretractionMaxVSelect,2);
SPEEDretractionMaxVNorm = normalize(SPEEDretractionMaxVSelect,2);
AccelretractionMaxVNorm = normalize(AccelretractionMaxVSelect,2);
FORCEretractionMaxVNorm = normalize(FORCEretractionMaxVSelect,2);

for i1 =  TimePoints
    protSegmentsNorm{i1} = normalize(protSegmentsSelect{i1},2);
    protSegmentsSpeedNorm{i1} = normalize(protSegmentsSpeedSelect{i1},2);
    protSegmentsAccelNorm{i1} = normalize(protSegmentsAccelSelect{i1},2);
    protSegmentsTfNorm{i1} = normalize(protSegmentsTfSelect{i1},2);
    retrSegmentsNorm{i1} = normalize(retrSegmentsSelect{i1},2);
    retrSegmentsSpeedNorm{i1} = normalize(retrSegmentsSpeedSelect{i1},2);
    retrSegmentsAccelNorm{i1} = normalize(retrSegmentsAccelSelect{i1},2);
    retrSegmentsTfNorm{i1} = normalize(retrSegmentsTfSelect{i1},2);
end



%% feather plot

nonPA = {InteDretractionMaxVSelect,SPEEDretractionMaxVSelect,AccelretractionMaxVSelect,FORCEretractionMaxVSelect,InteDprotrusionMaxVSelect,SPEEDprotrusionMaxVSelect,AccelprotrusionMaxVSelect,FORCEprotrusionMaxVSelect};
PA = {retrSegmentsSelect,retrSegmentsSpeedSelect,retrSegmentsAccelSelect,retrSegmentsTfSelect,protSegmentsSelect,protSegmentsSpeedSelect,protSegmentsAccelSelect,protSegmentsTfSelect};
PlotArray_ylabel1 = {'Retraction distance','Retraction velocity','Retraction acceleration','Traction force','Protrusion distance','Protrusion velocity','Protrusion acceleration','Traction force'};
PlotArray_ylabel2 = {'(\mum)','(\mum/min)',' ','({\itPa})','(\mum)','(\mum/min)',' ','({\itPa})'};
PlotArray_title = {'PA-Rac1(T17N) retraction distance','PA-Rac1(T17N) retraction velocity','PA-Rac1(T17N) retraction acceleration','PA-Rac1(T17N) retraction force','PA-Rac1(T17N) protrusion distance','PA-Rac1(T17N) protrusion velocity','PA-Rac1(T17N) protrusion acceleration','PA-Rac1(T17N) protrusion force'};
PlotArray_legendLoc = {'northwest','southwest','northwest','northwest','northwest','northwest','southwest','southwest'};
PlotArray_peakT = {5,5,5,5,7,7,7,7};
PlotArray_SelectedT = {[58,59],[58,59],[58,59],[58,59],[56,57,58,59],[56,57,58,59],[56,57,58,59],[56,57,58,59]};
PlotArray_T = 56:59;
PlotArray_color = {'r','c','b','m'};
PlotArray_UnitFactor = {0.2,2.4,1,1,0.2,2.4,1,1};
PlotArray_legendPos = {[0.3 0.68 0.1 0.1],[0.4 0.38 0.1 0.1],[0.32 0.68 0.1 0.1],[0.32 0.68 0.1 0.1],[0.3 0.6 0.1 0.1],[0.57 0.38 0.1 0.1],[0.32 0.68 0.1 0.1],[0.32 0.68 0.1 0.1]};



% PlotArray_range for TPSelected
PlotArray_range = {[-25 35 0 0.45],[-25 35 -0.2 0.75],[-25 35 0 0.2],[-25 95 -25,25],[-35 25 0 0.8],[-35 25 0 1.2],[-35 25 -0.15 0.05],[-35 85 -30,30]};


for i = 1 : size(PA,2)
    figure(1)
    set(gcf,'color','w');
    
    plot((0-PlotArray_peakT{i}:31-PlotArray_peakT{i})*timeinterval,PlotArray_UnitFactor{i}*(nanmean(nonPA{i})-nanmean(nonPA{i}(:,1))),'o-k','LineWidth',5,'MarkerSize',15)
    hold off
    for i2 = PlotArray_SelectedT{i}
        hold on
        plot((61-i2-PlotArray_peakT{i}:size(PlotRange,2)-PlotArray_peakT{i})*timeinterval,PlotArray_UnitFactor{i}*(nanmean(PA{i}{1, i2}(:, 61-i2+1:size(PlotRange,2)+1)) - nanmean(PA{i}{1, i2}(:, 61-i2+1)) +(nanmean(nonPA{i}(:,61-i2+1))-nanmean(nonPA{i}(:,1)))),['.-' PlotArray_color{i2-PlotArray_SelectedT{i}(1)+1}],'LineWidth',5,'MarkerSize',25);
        hold off
        label{i2-PlotArray_SelectedT{i}(1)+1} = ['PA-start at ' num2str((61-i2-PlotArray_peakT{i})*timeinterval) ' s'];
    end
    
    axis(PlotArray_range{i})
    legend({'control',label{1:size(PlotArray_SelectedT{i},2)}},'Position',PlotArray_legendPos{i},'FontSize',34);
    legend boxoff
    xlabel('t(s) aligned by V_m_a_x','FontWeight','bold','FontSize',15)
    ylabel({PlotArray_ylabel1{i};PlotArray_ylabel2{i}},'FontWeight','bold','FontSize',15)
    title(PlotArray_title{i},'FontSize',15,'FontWeight','bold');
    set(gca,'FontSize',45,'FontWeight','bold');
    figure(1);  set(gcf,'Position',[0 0 2000 1500]); export_fig(gcf,[filepath(InitialArray) savefolder PlotArray_title{i} '_feather plot_TPSelected.fig']);
    figure(1);  set(gcf,'Position',[0 0 2000 1500]); export_fig(gcf,[filepath(InitialArray) savefolder PlotArray_title{i} '_feather plot_TPSelected.png']);
    close(1)
end


    
%% false discovery rate (FDR) test

for i = 1 : size(PA,2)
    tempT1{i} = nonPA{i}-nonPA{i}(:,1);
    for i2 = PlotArray_SelectedT{i}
        clearvars dataT1 dataT2 lengthT1 lengthT2 dataT1T2 sdeALL sdeRandT1 sdeRandT2 sdeRandT1T2 sdeRandALL pvaluesRandT1T2 sdeT1 sdeT2 sdeT1T2 pvaluesT1T2
        tempT2{i}{i2} = PA{i}{1, i2}(:, 61-i2+1:size(PlotRange,2)+1)- PA{i}{1, i2}(:, 61-i2+1) +(nanmean(nonPA{i}(:,61-i2+1))-nanmean(nonPA{i}(:,1)));
        dataT1 = tempT1{i}(2:end,62-i2:13);
        dataT2 = tempT2{i}{i2}(:,1:13-(62-i2)+1);
        lengthT1 = size(dataT1,1);
        lengthT2 = size(dataT2,1);
        dataT1T2 = [dataT1;dataT2];
        lengthT1T2min = min([lengthT1,lengthT2]);
%         sdeALL = nansum(abs(dataT1T2-nanmean(dataT1T2)),2);
        for i3 = 1 : 10000
            tempSeq = randperm(lengthT1+lengthT2);
            sdeRandT1 = nansum((dataT1T2(tempSeq(1:lengthT1T2min),:)-nanmean(dataT1T2(tempSeq(1:lengthT1T2min),:))).^2,2);
            sdeRandT2 = nansum((dataT1T2(tempSeq(lengthT1T2min+1:lengthT1T2min*2),:)-nanmean(dataT1T2(tempSeq(lengthT1T2min+1:lengthT1T2min*2),:))).^2,2);
            sdeRandT1T2(i3,:) = [sdeRandT1;sdeRandT2];
            sdeRandALL(i3,:) = nansum((dataT1T2(tempSeq(1:lengthT1T2min*2),:)-nanmean(dataT1T2(tempSeq(1:lengthT1T2min*2),:))).^2,2);
%             [hvalues(i3), pvalues(i3)] = ttest2(sdeRandT1T2(i3,:),sdeALL);
            
            tempSeqT1 = randperm(lengthT1);
            tempSeqT2 = randperm(lengthT2);
            sdeT1 = nansum((dataT1(tempSeqT1(1:lengthT1T2min),:)-nanmean(dataT1(tempSeqT1(1:lengthT1T2min),:))).^2,2);
            sdeT2 = nansum((dataT2(tempSeqT2(1:lengthT1T2min),:)-nanmean(dataT2(tempSeqT2(1:lengthT1T2min),:))).^2,2);
            sdeT1T2(i3,:) = [sdeT1;sdeT2];
            sdeALL(i3,:) = nansum(([dataT1(tempSeqT1(1:lengthT1T2min),:);dataT2(tempSeqT2(1:lengthT1T2min),:)]-nanmean([dataT1(tempSeqT1(1:lengthT1T2min),:);dataT2(tempSeqT2(1:lengthT1T2min),:)])).^2,2);
        end
        pvaluesRandT1T2 = mattest(sdeRandT1T2,sdeRandALL,'permute',true);
%         [fdr_RandT1T2,q_RandT1T2,priori_RandT1T2,R2_RandT1T2] = mafdr(pvaluesRandT1T2,'Method','polynomial');
        [fdr_RandT1T2,q_RandT1T2] = mafdr(pvaluesRandT1T2);
        
%         sdeT1 = nansum(abs(dataT1-nanmean(dataT1)),2);
%         sdeT2 = nansum(abs(dataT2-nanmean(dataT2)),2); 
%         sdeT1T2 = [sdeT1;sdeT2];
        pvaluesT1T2 = mattest(sdeT1T2,sdeALL,'permute',true);
%         [fdr_T1T2,q_T1T2,priori_T1T2,R2_T1T2] = mafdr(pvaluesT1T2,'Method','polynomial');
        [fdr_T1T2,q_T1T2] = mafdr(pvaluesT1T2);
        temp = sort(q_RandT1T2);
        if nanmean(q_T1T2) < temp(floor(length(temp)))
            h_results{i,i2} = 'significant';
        else if nanmean(q_T1T2) > temp(floor(length(temp)))
                h_results{i,i2} = 'non-significant';
            end
        end
        Q(i,i2) = nanmean(q_T1T2);
        Qrand(i,i2) = temp(floor(length(temp)));
    end
    
end


%% save data
xlswrite([filepath(InitialArray) savefolder 'Qvalue.xls'],Q,1,'A1')
xlswrite([filepath(InitialArray) savefolder 'Qvalue.xls'],Qrand,1,'A10')


save([filepath(InitialArray) savefolder 'AllData.mat'])