clear
clc

tic

MSGID = 'SPLINES:CHCKXYWP:NaNs';
warning('off', MSGID)

timeinterval = 10;

filepathALL{1} = 'Y:\GTPase4\HT1080RhoA2Gin1%FBSonPAA20190820\analysis\Windowing\';
SelectedCell{1} = [2, 3, 4, 8];

filepathALL{2} = 'Y:\GTPase4\HT1080RhoA2Gin1%FBSonPAA20190619 - HresTF\analysis\Windowing\';
SelectedCell{2} = [1, 2, 3, 4, 5, 6];

filepathALL{3} = 'Y:\GTPase4\HT1080RhoA2Gin1%FBSonPAA20190920\analysis\Windowing\';
SelectedCell{3} = [1, 2, 3, 4];

filepathALL{4} = 'Y:\GTPase4\HT1080RhoA2Gin1%FBSonPAA20190923\analysis\Windowing\';
SelectedCell{4} = [1, 2, 3, 4, 5];

filepathALL{5} = 'Y:\GTPase4\HT1080RhoA2Gin1%FBSonPAA20191002\analysis\Windowing\';
SelectedCell{5} = [1, 2, 3, 4];

LongestMovie = [5, 3]; % point out the longest movie for initial matrix generation

filename = 'HT1080RhoA2G';

savepath = 'Y:\GTPase4\Combined Results\RhoA\AnalysisResults_withSmoothedProtrusion_AllDepth_v2_201216y\';
mkdir(savepath);

a = 6; % min time range
b = 5; % min value range
d = 1 : 10; % FRET&C2 depth
h = 1 : 10; % Force depth



clearvars -except MSGID timeinterval filepathALL SelectedCell LongestMovie filename savepath a b d h

% create initial matrix for time series storage
filepath = filepathALL{LongestMovie(1)};
load([filepath filename 'in1%fbsonpaa00' num2str(LongestMovie(2)) 'xy1\WindowingPackage\protrusion_samples\protrusion_samples.mat']);
n = size(protSamples.avgNormal,2);


timerange = n-12 : n+12;


% define parameter
distance = a;
valuerange = b;
FRETdepth = d;
FORCEdepth = h;
C2depth = d;




for i11 = FRETdepth
    FRETprotrusion{i11}(1,:) = nan(1,2*n);
    FRETretraction{i11}(1,:) = nan(1,2*n);
    FORCEprotrusion{i11}(1,:) = nan(1,2*n);
    FORCEretraction{i11}(1,:) = nan(1,2*n);
    C2protrusion{i11}(1,:) = nan(1,2*n);
    C2retraction{i11}(1,:) = nan(1,2*n);
    FRETprotrusionMaxV{i11}(1,:) = nan(1,2*n);
    FORCEprotrusionMaxV{i11}(1,:) = nan(1,2*n);
    FRETretractionMaxV{i11}(1,:) = nan(1,2*n);
    FORCEretractionMaxV{i11}(1,:) = nan(1,2*n);
end

InteDprotrusion(1,:) = nan(1,2*n);
InteDretraction(1,:) = nan(1,2*n);
InteDprotrusionMaxV(1,:) = nan(1,2*n);
InteDretractionMaxV(1,:) = nan(1,2*n);

SPEEDprotrusion(1,:) = nan(1,2*n);
SPEEDretraction(1,:) = nan(1,2*n);
SPEEDprotrusionMaxV(1,:) = nan(1,2*n);
SPEEDretractionMaxV(1,:) = nan(1,2*n);



for g = 1 : length(filepathALL)
    filepath = filepathALL{g};
    a1 = SelectedCell{g};
    
    for j = a1
        for k = 1 : 1
            load([filepath filename 'in1%fbsonpaa00' num2str(j) 'xy1\WindowingPackage\protrusion_samples\protrusion_samples.mat']);
            
            FRET = load([filepath filename 'in1%fbsonpaa00' num2str(j) 'xy1\WindowingPackage\window_sampling\Raw images - channel 1.mat']);
            FORCE = load([filepath filename 'in1%fbsonpaa00' num2str(j) 'xy1\WindowingPackage\window_sampling\Raw images - channel 2.mat']);
            C2 = load([filepath filename 'in1%fbsonpaa00' num2str(j) 'xy1\WindowingPackage\window_sampling\Raw images - channel 3.mat']);
            
            % Align the mismached timeframes by window sampling
            FORCEchosen = squeeze(FORCE.samples.avg(:,1,:));
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
                
                protSamples.avgNormal(:,i4) = circshift(protSamples.avgNormal(:,i4),shift{i4});
                for i6 = 1 : size(FRET.samples.avg,2)
                    FRET.samples.avg(:,i6,i4) = circshift(FRET.samples.avg(:,i6,i4),shift{i4});
                    FORCE.samples.avg(:,i6,i4) = circshift(FORCE.samples.avg(:,i6,i4),shift{i4});
                    C2.samples.avg(:,i6,i4) = circshift(C2.samples.avg(:,i6,i4),shift{i4});
                end
            end
            
            
            
            [m2,n2] = size(protSamples.avgNormal);
            
            for i0 = 1 : n2-1
                protIntegrated{g,j,k}(:,i0) = nansum(protSamples.avgNormal(:,1:i0),2);
            end
            
            
            clearvars pp p peakval peakloc peakwidth peakprominence valleyloc valleyval retractiondistanceval protrusiondistanceval peakvalfilted peaklocfilted
            
            % smooth FRET & FORCE & C2
            for i = 1 : length(protIntegrated{g,j,k}(~isnan(protIntegrated{g,j,k}(:,1)),1))
                for tempi = FRETdepth
                    tempN = find(~isnan(FRET.samples.avg(i,tempi,:)));
                    if length(tempN)>2
                        % smooth FRET
                        [temp1,temp2] = csaps(1:length(tempN),FRET.samples.avg(i,tempi,tempN),0.9);
                        FRET.samples.avg(i,tempi,tempN(1:end-1)) = temp1.coefs(:,4);
                        % smooth FORCE
                        [temp1,temp2] = csaps(1:length(tempN),FORCE.samples.avg(i,tempi,tempN),0.9);
                        FORCE.samples.avg(i,tempi,tempN(1:end-1)) = temp1.coefs(:,4);
                        % smooth C2
                        [temp1,temp2] = csaps(1:length(tempN),C2.samples.avg(i,tempi,tempN),0.9);
                        C2.samples.avg(i,tempi,tempN(1:end-1)) = temp1.coefs(:,4);
                    end
                end
            end
            
            
            for i = 1 : length(protIntegrated{g,j,k}(~isnan(protIntegrated{g,j,k}(:,1)),1))
                tempN = find(~isnan(protIntegrated{g,j,k}(i,:)));
                if length(tempN)>2
                    [pp{i},p(i)] = csaps(1:n2-1 , protIntegrated{g,j,k}(i,:),0.1);
                    % get smoothed membrane dynamic & speed
                    protIntegratedSmoothed{g,j,k}(i,tempN(1:end-1)) = pp{i}.coefs(:,4);
                    protSamplesSmoothed{g,j,k}(i,1) = pp{i}.coefs(1,4);
                    for k1 = 2 : length(tempN)-1
                        protSamplesSmoothed{g,j,k}(i,tempN(k1)) = pp{i}.coefs(k1,4) - pp{i}.coefs(k1-1,4);
                    end
                    if size(pp{i}.coefs,2) == 4
                        [peakval{i}, peakloc{i}, peakwidth{i}, peakprominence{i}] = findpeaks(pp{i}.coefs(:,4));
                        
                        
                        % find local valley between each peak pairs
                        for i1 = 1 : length(peakval{i})-1
                            [minval, minidx] = min(pp{i}.coefs(peakloc{i}(i1):peakloc{i}(i1+1),4));
                            valleyloc{i}(i1) = peakloc{i}(i1) + minidx - 1;
                            valleyval{i}(i1) = minval;
                            retractiondistanceval(i,i1) = peakval{i}(i1) - valleyval{i}(i1);
                            protrusiondistanceval(i,i1) = peakval{i}(i1+1) - valleyval{i}(i1);
                        end
                        % remove the short term switches
                        peakvalfilted{i} = peakval{i};
                        peaklocfilted{i} = peakloc{i};
                        for i2 = 1 : n2-distance
                            temploc = find(peakloc{i}>i2 & peakloc{i}<=i2+distance);
                            temptemploc = find(peakval{i}(temploc) ~= max(peakval{i}(temploc)));
                            peaklocfilted{i}(temploc(temptemploc)) = NaN;
                            peakvalfilted{i}(temploc(temptemploc)) = NaN;
                        end
                        % get the protrusion & retraction arrays
                        peakposition = find(~isnan(peaklocfilted{i}));
                        for i3 = 2 : length(peakposition)-2
                            % add retraction array if its time length is larger than given distance and its range is larger than given scale.
                            if peakposition(i3) > 1 && valleyloc{i}(peakposition(i3)) - valleyloc{i}(peakposition(i3-1)) >= distance
                                if peakval{i}(peakposition(i3)) - valleyval{i}(peakposition(i3)) >= valuerange
                                    
                                    
                                    InteDretraction = [InteDretraction; nan(1,2*n)];
                                    SPEEDretraction = [SPEEDretraction; nan(1,2*n)];
                                    lengthbeforeretraction = peakloc{i}(peakposition(i3))-valleyloc{i}(peakposition(i3-1));
                                    
                                    InteDretraction(end,n-lengthbeforeretraction+1:n) = pp{1,i}.coefs(valleyloc{i}(peakposition(i3-1)):peakloc{i}(peakposition(i3))-1,4);
                                    SPEEDretraction(end,n-lengthbeforeretraction+1:n) = protSamplesSmoothed{g,j,k}(i,valleyloc{i}(peakposition(i3-1)):peakloc{i}(peakposition(i3))-1);
                                    
                                    lengthafterretraction = valleyloc{i}(peakposition(i3+1))-peakloc{i}(peakposition(i3))+1;
                                    
                                    InteDretraction(end,n+1:n+lengthafterretraction) = pp{1,i}.coefs(peakloc{i}(peakposition(i3)):valleyloc{i}(peakposition(i3+1)),4);
                                    SPEEDretraction(end,n+1:n+lengthafterretraction) = protSamplesSmoothed{g,j,k}(i,peakloc{i}(peakposition(i3)):valleyloc{i}(peakposition(i3+1)));
                                    
                                    
                                    InteDretractionMaxV = [InteDretractionMaxV ; nan(1,2*n)];
                                    SPEEDretractionMaxV = [SPEEDretractionMaxV ; nan(1,2*n)];
                                    clearvars retractionV pp1 p1 m1 n1
                                    retractionV = protSamples.avgNormal(i,valleyloc{i}(peakposition(i3-1)):valleyloc{i}(peakposition(i3+1)));
                                    [pp1,p1] = csaps(1:length(retractionV),retractionV,0.1);
                                    if size(pp1.coefs,2) == 4
                                        [m1,n1] = find(pp1.coefs(:,4) == min(pp1.coefs(:,4)));
                                        InteDretractionMaxV(end,n-m1+2:n) = pp{1,i}.coefs(valleyloc{i}(peakposition(i3-1)):valleyloc{i}(peakposition(i3-1))+m1-2,4);
                                        InteDretractionMaxV(end,n+1:n+length(retractionV)-m1+1) = pp{1,i}.coefs(valleyloc{i}(peakposition(i3-1))+m1-1:valleyloc{i}(peakposition(i3+1)),4);
                                        SPEEDretractionMaxV(end,n-m1+2:n) = protSamplesSmoothed{g,j,k}(i,valleyloc{i}(peakposition(i3-1)):valleyloc{i}(peakposition(i3-1))+m1-2);
                                        SPEEDretractionMaxV(end,n+1:n+length(retractionV)-m1+1) = protSamplesSmoothed{g,j,k}(i,valleyloc{i}(peakposition(i3-1))+m1-1:valleyloc{i}(peakposition(i3+1)));
                                    end
                                    
                                    
                                    for i12 = FRETdepth
                                        FRETretraction{i12} = [FRETretraction{i12} ; nan(1,2*n)];
                                        FORCEretraction{i12} = [FORCEretraction{i12} ; nan(1,2*n)];
                                        C2retraction{i12} = [C2retraction{i12}; nan(1,2*n)];
                                        
                                        lengthbeforeretraction = peakloc{i}(peakposition(i3))-valleyloc{i}(peakposition(i3-1));
                                        FRETretraction{i12}(end,n-lengthbeforeretraction+1:n) = FRET.samples.avg(i,i12,valleyloc{i}(peakposition(i3-1)):peakloc{i}(peakposition(i3))-1);
                                        FORCEretraction{i12}(end,n-lengthbeforeretraction+1:n) = FORCE.samples.avg(i,i12,valleyloc{i}(peakposition(i3-1)):peakloc{i}(peakposition(i3))-1);
                                        C2retraction{i12}(end,n-lengthbeforeretraction+1:n) = C2.samples.avg(i,i12,valleyloc{i}(peakposition(i3-1)):peakloc{i}(peakposition(i3))-1);
                                        
                                        
                                        lengthafterretraction = valleyloc{i}(peakposition(i3+1))-peakloc{i}(peakposition(i3))+1;
                                        FRETretraction{i12}(end,n+1:n+lengthafterretraction) = FRET.samples.avg(i,i12,peakloc{i}(peakposition(i3)):valleyloc{i}(peakposition(i3+1)));
                                        FORCEretraction{i12}(end,n+1:n+lengthafterretraction) = FORCE.samples.avg(i,i12,peakloc{i}(peakposition(i3)):valleyloc{i}(peakposition(i3+1)));
                                        C2retraction{i12}(end,n+1:n+lengthafterretraction) = C2.samples.avg(i,i12,peakloc{i}(peakposition(i3)):valleyloc{i}(peakposition(i3+1)));
                                        
                                        
                                        FRETretractionMaxV{i12} = [FRETretractionMaxV{i12} ; nan(1,2*n)];
                                        FORCEretractionMaxV{i12} = [FORCEretractionMaxV{i12} ; nan(1,2*n)];
                                        
                                        clearvars retractionV pp1 p1 m1 n1
                                        retractionV = protSamples.avgNormal(i,valleyloc{i}(peakposition(i3-1)):valleyloc{i}(peakposition(i3+1)));
                                        [pp1,p1] = csaps(1:length(retractionV),retractionV,0.1);
                                        if size(pp1.coefs,2) == 4
                                            [m1,n1] = find(pp1.coefs(:,4) == min(pp1.coefs(:,4)));
                                            FRETretractionMaxV{i12}(end,n-m1+2:n) = FRET.samples.avg(i,i12,valleyloc{i}(peakposition(i3-1)):valleyloc{i}(peakposition(i3-1))+m1-2);
                                            FORCEretractionMaxV{i12}(end,n-m1+2:n) = FORCE.samples.avg(i,i12,valleyloc{i}(peakposition(i3-1)):valleyloc{i}(peakposition(i3-1))+m1-2);
                                            FRETretractionMaxV{i12}(end,n+1:n+length(retractionV)-m1+1) = FRET.samples.avg(i,i12,valleyloc{i}(peakposition(i3-1))+m1-1:valleyloc{i}(peakposition(i3+1)));
                                            FORCEretractionMaxV{i12}(end,n+1:n+length(retractionV)-m1+1) = FORCE.samples.avg(i,i12,valleyloc{i}(peakposition(i3-1))+m1-1:valleyloc{i}(peakposition(i3+1)));
                                            
                                        end
                                        
                                    end
                                    
                                end
                                % add protrusion array if its length is larger than given distance.
                                if peakposition(i3) > 1 && peakloc{i}(peakposition(i3)) - peakloc{i}(peakposition(i3-1)) >= distance
                                    if peakval{i}(peakposition(i3)) - valleyval{i}(peakposition(i3-1)) >= valuerange
                                        
                                        
                                        InteDprotrusion = [InteDprotrusion; nan(1,2*n)];
                                        SPEEDprotrusion = [SPEEDprotrusion; nan(1,2*n)];
                                        lengthbeforeprotrusion = valleyloc{i}(peakposition(i3))-peakloc{i}(peakposition(i3-1));
                                        
                                        InteDprotrusion(end,n-lengthbeforeprotrusion+1:n) = pp{1,i}.coefs(peakloc{i}(peakposition(i3-1)):valleyloc{i}(peakposition(i3))-1,4);
                                        SPEEDprotrusion(end,n-lengthbeforeprotrusion+1:n) = protSamplesSmoothed{g,j,k}(i,peakloc{i}(peakposition(i3-1)):valleyloc{i}(peakposition(i3))-1);
                                        
                                        lengthafterprotrusion = peakloc{i}(peakposition(i3+1))-valleyloc{i}(peakposition(i3))+1;
                                        
                                        InteDprotrusion(end,n+1:n+lengthafterprotrusion) = pp{1,i}.coefs(valleyloc{i}(peakposition(i3)):peakloc{i}(peakposition(i3+1)),4);
                                        SPEEDprotrusion(end,n+1:n+lengthafterprotrusion) = protSamplesSmoothed{g,j,k}(i,valleyloc{i}(peakposition(i3)):peakloc{i}(peakposition(i3+1)));
                                        
                                        
                                        InteDprotrusionMaxV = [InteDprotrusionMaxV ; nan(1,2*n)];
                                        SPEEDprotrusionMaxV = [SPEEDprotrusionMaxV ; nan(1,2*n)];
                                        clearvars protrusionV pp1 p1 m1 n1
                                        protrusionV = protSamples.avgNormal(i,valleyloc{i}(peakposition(i3-1)):valleyloc{i}(peakposition(i3+1)));
                                        [pp1,p1] = csaps(1:length(protrusionV),protrusionV,0.1);
                                        if size(pp1.coefs,2) == 4
                                            [m1,n1] = find(pp1.coefs(:,4) == max(pp1.coefs(:,4)));
                                            InteDprotrusionMaxV(end,n-m1+2:n) = pp{1,i}.coefs(valleyloc{i}(peakposition(i3-1)):valleyloc{i}(peakposition(i3-1))+m1-2,4);
                                            InteDprotrusionMaxV(end,n+1:n+length(protrusionV)-m1+1) = pp{1,i}.coefs(valleyloc{i}(peakposition(i3-1))+m1-1:valleyloc{i}(peakposition(i3+1)),4);
                                            SPEEDprotrusionMaxV(end,n-m1+2:n) = protSamplesSmoothed{g,j,k}(i,valleyloc{i}(peakposition(i3-1)):valleyloc{i}(peakposition(i3-1))+m1-2);
                                            SPEEDprotrusionMaxV(end,n+1:n+length(protrusionV)-m1+1) = protSamplesSmoothed{g,j,k}(i,valleyloc{i}(peakposition(i3-1))+m1-1:valleyloc{i}(peakposition(i3+1)));
                                            
                                        end
                                        
                                        
                                        for i13 = FRETdepth
                                            FRETprotrusion{i13} = [FRETprotrusion{i13} ; nan(1,2*n)];
                                            FORCEprotrusion{i13} = [FORCEprotrusion{i13} ; nan(1,2*n)];
                                            C2protrusion{i13} = [C2protrusion{i13}; nan(1,2*n)];
                                            
                                            lengthbeforeprotrusion = valleyloc{i}(peakposition(i3))-peakloc{i}(peakposition(i3-1));
                                            FRETprotrusion{i13}(end,n-lengthbeforeprotrusion+1:n) = FRET.samples.avg(i,i13,peakloc{i}(peakposition(i3-1)):valleyloc{i}(peakposition(i3))-1);
                                            FORCEprotrusion{i13}(end,n-lengthbeforeprotrusion+1:n) = FORCE.samples.avg(i,i13,peakloc{i}(peakposition(i3-1)):valleyloc{i}(peakposition(i3))-1);
                                            C2protrusion{i13}(end,n-lengthbeforeprotrusion+1:n) = C2.samples.avg(i,i13,peakloc{i}(peakposition(i3-1)):valleyloc{i}(peakposition(i3))-1);
                                            
                                            
                                            lengthafterprotrusion = peakloc{i}(peakposition(i3+1))-valleyloc{i}(peakposition(i3))+1;
                                            FRETprotrusion{i13}(end,n+1:n+lengthafterprotrusion) = FRET.samples.avg(i,i13,valleyloc{i}(peakposition(i3)):peakloc{i}(peakposition(i3+1)));
                                            FORCEprotrusion{i13}(end,n+1:n+lengthafterprotrusion) = FORCE.samples.avg(i,i13,valleyloc{i}(peakposition(i3)):peakloc{i}(peakposition(i3+1)));
                                            C2protrusion{i13}(end,n+1:n+lengthafterprotrusion) = C2.samples.avg(i,i13,valleyloc{i}(peakposition(i3)):peakloc{i}(peakposition(i3+1)));
                                            
                                            
                                            FRETprotrusionMaxV{i13} = [FRETprotrusionMaxV{i13} ; nan(1,2*n)];
                                            FORCEprotrusionMaxV{i13} = [FORCEprotrusionMaxV{i13} ; nan(1,2*n)];
                                            
                                            clearvars protrusionV pp1 p1 m1 n1
                                            protrusionV = protSamples.avgNormal(i,valleyloc{i}(peakposition(i3-1)):valleyloc{i}(peakposition(i3+1)));
                                            [pp1,p1] = csaps(1:length(protrusionV),protrusionV,0.1);
                                            if size(pp1.coefs,2) == 4
                                                [m1,n1] = find(pp1.coefs(:,4) == max(pp1.coefs(:,4)));
                                                FRETprotrusionMaxV{i13}(end,n-m1+2:n) = FRET.samples.avg(i,i13,valleyloc{i}(peakposition(i3-1)):valleyloc{i}(peakposition(i3-1))+m1-2);
                                                FORCEprotrusionMaxV{i13}(end,n-m1+2:n) = FORCE.samples.avg(i,i13,valleyloc{i}(peakposition(i3-1)):valleyloc{i}(peakposition(i3-1))+m1-2);
                                                FRETprotrusionMaxV{i13}(end,n+1:n+length(protrusionV)-m1+1) = FRET.samples.avg(i,i13,valleyloc{i}(peakposition(i3-1))+m1-1:valleyloc{i}(peakposition(i3+1)));
                                                FORCEprotrusionMaxV{i13}(end,n+1:n+length(protrusionV)-m1+1) = FORCE.samples.avg(i,i13,valleyloc{i}(peakposition(i3-1))+m1-1:valleyloc{i}(peakposition(i3+1)));
                                                
                                            end
                                            
                                        end
                                    end
                                end
                            end
                            %                             [j k i toc]
                        end
                        
                        
                        % cell based data normalization
                        [ProM((j-1)*6+k),ProN((j-1)*6+k)] = size(FRETprotrusion{1});
                        [RetM((j-1)*6+k),RetN((j-1)*6+k)] = size(FRETretraction{1});
                        
                        for i14 = FRETdepth
                            if j == 1 && k == 1
                                X = FRETretraction{i14}(1:RetM((j-1)*6+k),:);
                                xmu=nanmean(X);
                                xsigma=nanstd(X);
                                FRETretractionCellNorm{i14}(1:RetM((j-1)*6+k),:) = (X-repmat(xmu,size(X,1),1))./repmat(xsigma,size(X,1),1);
                                X = FORCEretraction{i14}(1:RetM((j-1)*6+k),:);
                                xmu=nanmean(X);
                                xsigma=nanstd(X);
                                FORCEretractionCellNorm{i14}(1:RetM((j-1)*6+k),:) = (X-repmat(xmu,size(X,1),1))./repmat(xsigma,size(X,1),1);
                                X = C2retraction{i14}(1:RetM((j-1)*6+k),:);
                                xmu=nanmean(X);
                                xsigma=nanstd(X);
                                C2retractionCellNorm{i14}(1:RetM((j-1)*6+k),:) = (X-repmat(xmu,size(X,1),1))./repmat(xsigma,size(X,1),1);
                                X = FRETprotrusion{i14}(1:ProM((j-1)*6+k),:);
                                xmu=nanmean(X);
                                xsigma=nanstd(X);
                                FRETprotrusionCellNorm{i14}(1:ProM((j-1)*6+k),:) = (X-repmat(xmu,size(X,1),1))./repmat(xsigma,size(X,1),1);
                                X = FORCEprotrusion{i14}(1:ProM((j-1)*6+k),:);
                                xmu=nanmean(X);
                                xsigma=nanstd(X);
                                FORCEprotrusionCellNorm{i14}(1:ProM((j-1)*6+k),:) = (X-repmat(xmu,size(X,1),1))./repmat(xsigma,size(X,1),1);
                                X = C2protrusion{i14}(1:ProM((j-1)*6+k),:);
                                xmu=nanmean(X);
                                xsigma=nanstd(X);
                                C2protrusionCellNorm{i14}(1:ProM((j-1)*6+k),:) = (X-repmat(xmu,size(X,1),1))./repmat(xsigma,size(X,1),1);
                                X = FRETretractionMaxV{i14}(1:RetM((j-1)*6+k),:);
                                xmu=nanmean(X);
                                xsigma=nanstd(X);
                                FRETretractionMaxVCellNorm{i14}(1:RetM((j-1)*6+k),:) = (X-repmat(xmu,size(X,1),1))./repmat(xsigma,size(X,1),1);
                                X = FORCEretractionMaxV{i14}(1:RetM((j-1)*6+k),:);
                                xmu=nanmean(X);
                                xsigma=nanstd(X);
                                FORCEretractionMaxVCellNorm{i14}(1:RetM((j-1)*6+k),:) = (X-repmat(xmu,size(X,1),1))./repmat(xsigma,size(X,1),1);
                                X = FRETprotrusionMaxV{i14}(1:ProM((j-1)*6+k),:);
                                xmu=nanmean(X);
                                xsigma=nanstd(X);
                                FRETprotrusionMaxVCellNorm{i14}(1:ProM((j-1)*6+k),:) = (X-repmat(xmu,size(X,1),1))./repmat(xsigma,size(X,1),1);
                                X = FORCEprotrusionMaxV{i14}(1:ProM((j-1)*6+k),:);
                                xmu=nanmean(X);
                                xsigma=nanstd(X);
                                FORCEprotrusionMaxVCellNorm{i14}(1:ProM((j-1)*6+k),:) = (X-repmat(xmu,size(X,1),1))./repmat(xsigma,size(X,1),1);
                                X = InteDprotrusionMaxV(1:ProM((j-1)*6+k),:);
                                xmu=nanmean(X);
                                xsigma=nanstd(X);
                                InteDprotrusionMaxVCellNorm(1:ProM((j-1)*6+k),:) = (X-repmat(xmu,size(X,1),1))./repmat(xsigma,size(X,1),1);
                                X = InteDprotrusion(1:ProM((j-1)*6+k),:);
                                xmu=nanmean(X);
                                xsigma=nanstd(X);
                                InteDprotrusionCellNorm(1:ProM((j-1)*6+k),:) = (X-repmat(xmu,size(X,1),1))./repmat(xsigma,size(X,1),1);
                                X = InteDretractionMaxV(1:RetM((j-1)*6+k),:);
                                xmu=nanmean(X);
                                xsigma=nanstd(X);
                                InteDretractionMaxVCellNorm(1:RetM((j-1)*6+k),:) = (X-repmat(xmu,size(X,1),1))./repmat(xsigma,size(X,1),1);
                                X = InteDretraction(1:RetM((j-1)*6+k),:);
                                xmu=nanmean(X);
                                xsigma=nanstd(X);
                                InteDretractionCellNorm(1:RetM((j-1)*6+k),:) = (X-repmat(xmu,size(X,1),1))./repmat(xsigma,size(X,1),1);
                                X = SPEEDprotrusionMaxV(1:ProM((j-1)*6+k),:);
                                xmu=nanmean(X);
                                xsigma=nanstd(X);
                                SPEEDprotrusionMaxVCellNorm(1:ProM((j-1)*6+k),:) = (X-repmat(xmu,size(X,1),1))./repmat(xsigma,size(X,1),1);
                                X = SPEEDprotrusion(1:ProM((j-1)*6+k),:);
                                xmu=nanmean(X);
                                xsigma=nanstd(X);
                                SPEEDprotrusionCellNorm(1:ProM((j-1)*6+k),:) = (X-repmat(xmu,size(X,1),1))./repmat(xsigma,size(X,1),1);
                                X = SPEEDretractionMaxV(1:RetM((j-1)*6+k),:);
                                xmu=nanmean(X);
                                xsigma=nanstd(X);
                                SPEEDretractionMaxVCellNorm(1:RetM((j-1)*6+k),:) = (X-repmat(xmu,size(X,1),1))./repmat(xsigma,size(X,1),1);
                                X = SPEEDretraction(1:RetM((j-1)*6+k),:);
                                xmu=nanmean(X);
                                xsigma=nanstd(X);
                                SPEEDretractionCellNorm(1:RetM((j-1)*6+k),:) = (X-repmat(xmu,size(X,1),1))./repmat(xsigma,size(X,1),1);
                            else
                                X = FRETretraction{i14}(RetM((j-1)*6+k-1)+1:RetM((j-1)*6+k),:);
                                xmu=nanmean(X);
                                xsigma=nanstd(X);
                                FRETretractionCellNorm{i14}(RetM((j-1)*6+k-1)+1:RetM((j-1)*6+k),:) = (X-repmat(xmu,size(X,1),1))./repmat(xsigma,size(X,1),1);
                                X = FORCEretraction{i14}(RetM((j-1)*6+k-1)+1:RetM((j-1)*6+k),:);
                                xmu=nanmean(X);
                                xsigma=nanstd(X);
                                FORCEretractionCellNorm{i14}(RetM((j-1)*6+k-1)+1:RetM((j-1)*6+k),:) = (X-repmat(xmu,size(X,1),1))./repmat(xsigma,size(X,1),1);
                                X = C2retraction{i14}(RetM((j-1)*6+k-1)+1:RetM((j-1)*6+k),:);
                                xmu=nanmean(X);
                                xsigma=nanstd(X);
                                C2retractionCellNorm{i14}(RetM((j-1)*6+k-1)+1:RetM((j-1)*6+k),:) = (X-repmat(xmu,size(X,1),1))./repmat(xsigma,size(X,1),1);
                                X = FRETprotrusion{i14}(ProM((j-1)*6+k-1)+1:ProM((j-1)*6+k),:);
                                xmu=nanmean(X);
                                xsigma=nanstd(X);
                                FRETprotrusionCellNorm{i14}(ProM((j-1)*6+k-1)+1:ProM((j-1)*6+k),:) = (X-repmat(xmu,size(X,1),1))./repmat(xsigma,size(X,1),1);
                                X = FORCEprotrusion{i14}(ProM((j-1)*6+k-1)+1:ProM((j-1)*6+k),:);
                                xmu=nanmean(X);
                                xsigma=nanstd(X);
                                FORCEprotrusionCellNorm{i14}(ProM((j-1)*6+k-1)+1:ProM((j-1)*6+k),:) = (X-repmat(xmu,size(X,1),1))./repmat(xsigma,size(X,1),1);
                                X = C2protrusion{i14}(ProM((j-1)*6+k-1)+1:ProM((j-1)*6+k),:);
                                xmu=nanmean(X);
                                xsigma=nanstd(X);
                                C2protrusionCellNorm{i14}(ProM((j-1)*6+k-1)+1:ProM((j-1)*6+k),:) = (X-repmat(xmu,size(X,1),1))./repmat(xsigma,size(X,1),1);
                                X = FRETretractionMaxV{i14}(RetM((j-1)*6+k-1)+1:RetM((j-1)*6+k),:);
                                xmu=nanmean(X);
                                xsigma=nanstd(X);
                                FRETretractionMaxVCellNorm{i14}(RetM((j-1)*6+k-1)+1:RetM((j-1)*6+k),:) = (X-repmat(xmu,size(X,1),1))./repmat(xsigma,size(X,1),1);
                                X = FORCEretractionMaxV{i14}(RetM((j-1)*6+k-1)+1:RetM((j-1)*6+k),:);
                                xmu=nanmean(X);
                                xsigma=nanstd(X);
                                FORCEretractionMaxVCellNorm{i14}(RetM((j-1)*6+k-1)+1:RetM((j-1)*6+k),:) = (X-repmat(xmu,size(X,1),1))./repmat(xsigma,size(X,1),1);
                                X = FRETprotrusionMaxV{i14}(ProM((j-1)*6+k-1)+1:ProM((j-1)*6+k),:);
                                xmu=nanmean(X);
                                xsigma=nanstd(X);
                                FRETprotrusionMaxVCellNorm{i14}(ProM((j-1)*6+k-1)+1:ProM((j-1)*6+k),:) = (X-repmat(xmu,size(X,1),1))./repmat(xsigma,size(X,1),1);
                                X = FORCEprotrusionMaxV{i14}(ProM((j-1)*6+k-1)+1:ProM((j-1)*6+k),:);
                                xmu=nanmean(X);
                                xsigma=nanstd(X);
                                FORCEprotrusionMaxVCellNorm{i14}(ProM((j-1)*6+k-1)+1:ProM((j-1)*6+k),:) = (X-repmat(xmu,size(X,1),1))./repmat(xsigma,size(X,1),1);
                                X = InteDprotrusionMaxV(ProM((j-1)*6+k-1)+1:ProM((j-1)*6+k),:);
                                xmu=nanmean(X);
                                xsigma=nanstd(X);
                                InteDprotrusionMaxVCellNorm(ProM((j-1)*6+k-1)+1:ProM((j-1)*6+k),:) = (X-repmat(xmu,size(X,1),1))./repmat(xsigma,size(X,1),1);
                                X = InteDprotrusion(ProM((j-1)*6+k-1)+1:ProM((j-1)*6+k),:);
                                xmu=nanmean(X);
                                xsigma=nanstd(X);
                                InteDprotrusionCellNorm(ProM((j-1)*6+k-1)+1:ProM((j-1)*6+k),:) = (X-repmat(xmu,size(X,1),1))./repmat(xsigma,size(X,1),1);
                                X = InteDretractionMaxV(RetM((j-1)*6+k-1)+1:RetM((j-1)*6+k),:);
                                xmu=nanmean(X);
                                xsigma=nanstd(X);
                                InteDretractionMaxVCellNorm(RetM((j-1)*6+k-1)+1:RetM((j-1)*6+k),:) = (X-repmat(xmu,size(X,1),1))./repmat(xsigma,size(X,1),1);
                                X = InteDretraction(RetM((j-1)*6+k-1)+1:RetM((j-1)*6+k),:);
                                xmu=nanmean(X);
                                xsigma=nanstd(X);
                                InteDretractionCellNorm(RetM((j-1)*6+k-1)+1:RetM((j-1)*6+k),:) = (X-repmat(xmu,size(X,1),1))./repmat(xsigma,size(X,1),1);
                                X = SPEEDprotrusionMaxV(ProM((j-1)*6+k-1)+1:ProM((j-1)*6+k),:);
                                xmu=nanmean(X);
                                xsigma=nanstd(X);
                                SPEEDprotrusionMaxVCellNorm(ProM((j-1)*6+k-1)+1:ProM((j-1)*6+k),:) = (X-repmat(xmu,size(X,1),1))./repmat(xsigma,size(X,1),1);
                                X = SPEEDprotrusion(ProM((j-1)*6+k-1)+1:ProM((j-1)*6+k),:);
                                xmu=nanmean(X);
                                xsigma=nanstd(X);
                                SPEEDprotrusionCellNorm(ProM((j-1)*6+k-1)+1:ProM((j-1)*6+k),:) = (X-repmat(xmu,size(X,1),1))./repmat(xsigma,size(X,1),1);
                                X = SPEEDretractionMaxV(RetM((j-1)*6+k-1)+1:RetM((j-1)*6+k),:);
                                xmu=nanmean(X);
                                xsigma=nanstd(X);
                                SPEEDretractionMaxVCellNorm(RetM((j-1)*6+k-1)+1:RetM((j-1)*6+k),:) = (X-repmat(xmu,size(X,1),1))./repmat(xsigma,size(X,1),1);
                                X = SPEEDretraction(RetM((j-1)*6+k-1)+1:RetM((j-1)*6+k),:);
                                xmu=nanmean(X);
                                xsigma=nanstd(X);
                                SPEEDretractionCellNorm(RetM((j-1)*6+k-1)+1:RetM((j-1)*6+k),:) = (X-repmat(xmu,size(X,1),1))./repmat(xsigma,size(X,1),1);
                                
                            end
                            
                        end
                        
                    end
                end
            end
        end
        
        
    end
end




% line based data normalization
for i15 = FRETdepth
    FRETretractionLineNorm{i15} = normalize(FRETretraction{i15},2);
    FORCEretractionLineNorm{i15} = normalize(FORCEretraction{i15},2);
    FRETprotrusionLineNorm{i15} = normalize(FRETprotrusion{i15},2);
    FORCEprotrusionLineNorm{i15} = normalize(FORCEprotrusion{i15},2);
    C2protrusionLineNorm{i15} = normalize(C2protrusion{i15},2);
    C2retractionLineNorm{i15} = normalize(C2retraction{i15},2);
    FRETretractionMaxVLineNorm{i15} = normalize(FRETretractionMaxV{i15},2);
    FORCEretractionMaxVLineNorm{i15} = normalize(FORCEretractionMaxV{i15},2);
    FRETprotrusionMaxVLineNorm{i15} = normalize(FRETprotrusionMaxV{i15},2);
    FORCEprotrusionMaxVLineNorm{i15} = normalize(FORCEprotrusionMaxV{i15},2);
end
InteDretractionLineNorm = normalize(InteDretraction,2);
InteDretractionMaxVLineNorm = normalize(InteDretractionMaxV,2);
InteDprotrusionLineNorm = normalize(InteDprotrusion,2);
InteDprotrusionMaxVLineNorm = normalize(InteDprotrusionMaxV,2);
SPEEDretractionLineNorm = normalize(SPEEDretraction,2);
SPEEDretractionMaxVLineNorm = normalize(SPEEDretractionMaxV,2);
SPEEDprotrusionLineNorm = normalize(SPEEDprotrusion,2);
SPEEDprotrusionMaxVLineNorm = normalize(SPEEDprotrusionMaxV,2);


save([savepath 'Distance' num2str(a) 'Valuerange' num2str(b) 'FRETdepth1to15FORCEdepth1to15.mat']);

load([savepath 'Distance' num2str(a) 'Valuerange' num2str(b) 'FRETdepth1to15FORCEdepth1to15.mat']);

for i16 = FRETdepth
    FORCEprotrusionMean(i16,:) = nanmean(FORCEprotrusion{i16}(:,timerange+1));
    FORCEprotrusionMaxVMean(i16,:) = nanmean(FORCEprotrusionMaxV{i16}(:,timerange+1));
    FORCEretractionMean(i16,:) = nanmean(FORCEretraction{i16}(:,timerange+1));
    FORCEretractionMaxVMean(i16,:) = nanmean(FORCEretractionMaxV{i16}(:,timerange+1));
    
    FRETprotrusionMean(i16,:) = nanmean(FRETprotrusion{i16}(:,timerange+1));
    FRETprotrusionMaxVMean(i16,:) = nanmean(FRETprotrusionMaxV{i16}(:,timerange+1));
    FRETretractionMean(i16,:) = nanmean(FRETretraction{i16}(:,timerange+1));
    FRETretractionMaxVMean(i16,:) = nanmean(FRETretractionMaxV{i16}(:,timerange+1));
    
    FORCEprotrusionLineNormMean(i16,:) = nanmean(FORCEprotrusionLineNorm{i16}(:,timerange+1));
    FORCEprotrusionMaxVLineNormMean(i16,:) = nanmean(FORCEprotrusionMaxVLineNorm{i16}(:,timerange+1));
    FORCEretractionLineNormMean(i16,:) = nanmean(FORCEretractionLineNorm{i16}(:,timerange+1));
    FORCEretractionMaxVLineNormMean(i16,:) = nanmean(FORCEretractionMaxVLineNorm{i16}(:,timerange+1));
    
    FRETprotrusionLineNormMean(i16,:) = nanmean(FRETprotrusionLineNorm{i16}(:,timerange+1));
    FRETprotrusionMaxVLineNormMean(i16,:) = nanmean(FRETprotrusionMaxVLineNorm{i16}(:,timerange+1));
    FRETretractionLineNormMean(i16,:) = nanmean(FRETretractionLineNorm{i16}(:,timerange+1));
    FRETretractionMaxVLineNormMean(i16,:) = nanmean(FRETretractionMaxVLineNorm{i16}(:,timerange+1));
end




FRETprotrusionMeanNorm = (FRETprotrusionMean - min(FRETprotrusionMean,[],'all'))/(max(FRETprotrusionMean,[],'all') - min(FRETprotrusionMean,[],'all'));
FRETprotrusionMaxVMeanNorm = (FRETprotrusionMaxVMean - min(FRETprotrusionMaxVMean,[],'all'))/(max(FRETprotrusionMaxVMean,[],'all') - min(FRETprotrusionMaxVMean,[],'all'));
FRETretractionMeanNorm = (FRETretractionMean - min(FRETretractionMean,[],'all'))/(max(FRETretractionMean,[],'all') - min(FRETretractionMean,[],'all'));
FRETretractionMaxVMeanNorm = (FRETretractionMaxVMean - min(FRETretractionMaxVMean,[],'all'))/(max(FRETretractionMaxVMean,[],'all') - min(FRETretractionMaxVMean,[],'all'));

FORCEprotrusionMeanNorm = (FORCEprotrusionMean - min(FORCEprotrusionMean,[],'all'))/(max(FORCEprotrusionMean,[],'all') - min(FORCEprotrusionMean,[],'all'));
FORCEprotrusionMaxVMeanNorm = (FORCEprotrusionMaxVMean - min(FORCEprotrusionMaxVMean,[],'all'))/(max(FORCEprotrusionMaxVMean,[],'all') - min(FORCEprotrusionMaxVMean,[],'all'));
FORCEretractionMeanNorm = (FORCEretractionMean - min(FORCEretractionMean,[],'all'))/(max(FORCEretractionMean,[],'all') - min(FORCEretractionMean,[],'all'));
FORCEretractionMaxVMeanNorm = (FORCEretractionMaxVMean - min(FORCEretractionMaxVMean,[],'all'))/(max(FORCEretractionMaxVMean,[],'all') - min(FORCEretractionMaxVMean,[],'all'));

FRETprotrusionLineNormMeanNorm = (FRETprotrusionLineNormMean - min(FRETprotrusionLineNormMean,[],'all'))/(max(FRETprotrusionLineNormMean,[],'all') - min(FRETprotrusionLineNormMean,[],'all'));
FRETprotrusionMaxVLineNormMeanNorm = (FRETprotrusionMaxVLineNormMean - min(FRETprotrusionMaxVLineNormMean,[],'all'))/(max(FRETprotrusionMaxVLineNormMean,[],'all') - min(FRETprotrusionMaxVLineNormMean,[],'all'));
FRETretractionLineNormMeanNorm = (FRETretractionLineNormMean - min(FRETretractionLineNormMean,[],'all'))/(max(FRETretractionLineNormMean,[],'all') - min(FRETretractionLineNormMean,[],'all'));
FRETretractionMaxVLineNormMeanNorm = (FRETretractionMaxVLineNormMean - min(FRETretractionMaxVLineNormMean,[],'all'))/(max(FRETretractionMaxVLineNormMean,[],'all') - min(FRETretractionMaxVLineNormMean,[],'all'));

FORCEprotrusionLineNormMeanNorm = (FORCEprotrusionLineNormMean - min(FORCEprotrusionLineNormMean,[],'all'))/(max(FORCEprotrusionLineNormMean,[],'all') - min(FORCEprotrusionLineNormMean,[],'all'));
FORCEprotrusionMaxVLineNormMeanNorm = (FORCEprotrusionMaxVLineNormMean - min(FORCEprotrusionMaxVLineNormMean,[],'all'))/(max(FORCEprotrusionMaxVLineNormMean,[],'all') - min(FORCEprotrusionMaxVLineNormMean,[],'all'));
FORCEretractionLineNormMeanNorm = (FORCEretractionLineNormMean - min(FORCEretractionLineNormMean,[],'all'))/(max(FORCEretractionLineNormMean,[],'all') - min(FORCEretractionLineNormMean,[],'all'));
FORCEretractionMaxVLineNormMeanNorm = (FORCEretractionMaxVLineNormMean - min(FORCEretractionMaxVLineNormMean,[],'all'))/(max(FORCEretractionMaxVLineNormMean,[],'all') - min(FORCEretractionMaxVLineNormMean,[],'all'));


halfT = ceil(length(timerange)/2);


for i17 = timerange-n
    figure(1)
    hold on;
    for i18 = FRETdepth
        if ~isnan(FORCEprotrusionLineNormMeanNorm(i18,i17+halfT))
            viscircles([0,0],length(FRETdepth)-i18+1,'color',[FORCEprotrusionLineNormMeanNorm(i18,i17+halfT),0,0],'LineWidth',30);
        end
    end
    viscircles([0,0],length(FRETdepth)+1,'color',[1,1,1],'LineWidth',35);
    hold off;
    axis equal
    axis([0 length(FRETdepth)+1 0 length(FRETdepth)+1]);
    title(['FORCEprotrusionLineNormMean' ' t=' num2str(i17)],'FontSize',20);
    annotation('textbox',[.55,.8,.1,.1],'string',['FORCE protrusion t=' num2str(i17)],'FontSize',15,'FontWeight','bold');
    windowMax = gcf; windowMax.WindowState = 'Maximize';
    pause(.1);
    M(i17+halfT) = getframe();
    close all
end
v = VideoWriter([savepath 'FORCEprotrusionLineNormMeanNorm.avi']);
v.FrameRate = 2;
open(v)
writeVideo(v,M)
close(v)

for i17 = timerange-n
    figure(1)
    hold on;
    for i18 = FRETdepth
        if ~isnan(FORCEprotrusionMaxVLineNormMeanNorm(i18,i17+halfT))
            viscircles([0,0],length(FRETdepth)-i18+1,'color',[FORCEprotrusionMaxVLineNormMeanNorm(i18,i17+halfT),0,0],'LineWidth',30);
        end
    end
    viscircles([0,0],length(FRETdepth)+1,'color',[1,1,1],'LineWidth',35);
    hold off;
    axis equal
    axis([0 length(FRETdepth)+1 0 length(FRETdepth)+1]);
    title(['FORCEprotrusionMaxVLineNormMean' ' t=' num2str(i17)],'FontSize',20);
    annotation('textbox',[.55,.8,.1,.1],'string',['FORCE protrusionMaxV t=' num2str(i17)],'FontSize',15,'FontWeight','bold');
    windowMax = gcf; windowMax.WindowState = 'Maximize';
    pause(.1);
    M(i17+halfT) = getframe();
    close all
end
v = VideoWriter([savepath 'FORCEprotrusionMaxVLineNormMeanNorm.avi']);
v.FrameRate = 2;
open(v)
writeVideo(v,M)
close(v)


for i17 = timerange-n
    figure(1)
    hold on;
    for i18 = FRETdepth
        if ~isnan(FRETprotrusionLineNormMeanNorm(i18,i17+halfT))
            viscircles([0,0],length(FRETdepth)-i18+1,'color',[0,FRETprotrusionLineNormMeanNorm(i18,i17+halfT),0],'LineWidth',30);
        end
    end
    viscircles([0,0],length(FRETdepth)+1,'color',[1,1,1],'LineWidth',35);
    hold off;
    axis equal
    axis([0 length(FRETdepth)+1 0 length(FRETdepth)+1]);
    title(['FRETprotrusionLineNormMean' ' t=' num2str(i17)],'FontSize',20);
    annotation('textbox',[.55,.8,.1,.1],'string',['FRET protrusion t=' num2str(i17)],'FontSize',15,'FontWeight','bold');
    windowMax = gcf; windowMax.WindowState = 'Maximize';
    pause(.1);
    M(i17+halfT) = getframe();
    close all
end
v = VideoWriter([savepath 'FRETprotrusionLineNormMeanNorm.avi']);
v.FrameRate = 2;
open(v)
writeVideo(v,M)
close(v)

for i17 = timerange-n
    figure(1)
    hold on;
    for i18 = FRETdepth
        if ~isnan(FRETprotrusionMaxVLineNormMeanNorm(i18,i17+halfT))
            viscircles([0,0],length(FRETdepth)-i18+1,'color',[0,FRETprotrusionMaxVLineNormMeanNorm(i18,i17+halfT),0],'LineWidth',30);
        end
    end
    viscircles([0,0],length(FRETdepth)+1,'color',[1,1,1],'LineWidth',35);
    hold off;
    axis equal
    axis([0 length(FRETdepth)+1 0 length(FRETdepth)+1]);
    title(['FRETprotrusionMaxVLineNormMean' ' t=' num2str(i17)],'FontSize',20);
    annotation('textbox',[.55,.8,.1,.1],'string',['FRET protrusionMaxV t=' num2str(i17)],'FontSize',15,'FontWeight','bold');
    windowMax = gcf; windowMax.WindowState = 'Maximize';
    pause(.1);
    M(i17+halfT) = getframe();
    close all
end
v = VideoWriter([savepath 'FRETprotrusionMaxVLineNormMeanNorm.avi']);
v.FrameRate = 2;
open(v)
writeVideo(v,M)
close(v)



for i17 = timerange-n
    figure(1)
    hold on;
    for i18 = FRETdepth
        if ~isnan(FORCEretractionLineNormMeanNorm(i18,i17+halfT))
            viscircles([0,0],length(FRETdepth)-i18+1,'color',[FORCEretractionLineNormMeanNorm(i18,i17+halfT),0,0],'LineWidth',30);
        end
    end
    viscircles([0,0],length(FRETdepth)+1,'color',[1,1,1],'LineWidth',35);
    hold off;
    axis equal
    axis([0 length(FRETdepth)+1 0 length(FRETdepth)+1]);
    title(['FORCEretractionLineNormMean' ' t=' num2str(i17)],'FontSize',20);
    annotation('textbox',[.55,.8,.1,.1],'string',['FORCE retraction t=' num2str(i17)],'FontSize',15,'FontWeight','bold');
    windowMax = gcf; windowMax.WindowState = 'Maximize';
    pause(.1);
    M(i17+halfT) = getframe();
    close all
end
v = VideoWriter([savepath 'FORCEretractionLineNormMeanNorm.avi']);
v.FrameRate = 2;
open(v)
writeVideo(v,M)
close(v)

for i17 = timerange-n
    figure(1)
    hold on;
    for i18 = FRETdepth
        if ~isnan(FORCEretractionMaxVLineNormMeanNorm(i18,i17+halfT))
            viscircles([0,0],length(FRETdepth)-i18+1,'color',[FORCEretractionMaxVLineNormMeanNorm(i18,i17+halfT),0,0],'LineWidth',30);
        end
    end
    viscircles([0,0],length(FRETdepth)+1,'color',[1,1,1],'LineWidth',35);
    hold off;
    axis equal
    axis([0 length(FRETdepth)+1 0 length(FRETdepth)+1]);
    title(['FORCEretractionMaxVLineNormMean' ' t=' num2str(i17)],'FontSize',20);
    annotation('textbox',[.55,.8,.1,.1],'string',['FORCE retractionMaxV t=' num2str(i17)],'FontSize',15,'FontWeight','bold');
    windowMax = gcf; windowMax.WindowState = 'Maximize';
    pause(.1);
    M(i17+halfT) = getframe();
    close all
end
v = VideoWriter([savepath 'FORCEretractionMaxVLineNormMeanNorm.avi']);
v.FrameRate = 2;
open(v)
writeVideo(v,M)
close(v)


for i17 = timerange-n
    figure(1)
    hold on;
    for i18 = FRETdepth
        if ~isnan(FRETretractionLineNormMeanNorm(i18,i17+halfT))
            viscircles([0,0],length(FRETdepth)-i18+1,'color',[0,FRETretractionLineNormMeanNorm(i18,i17+halfT),0],'LineWidth',30);
        end
    end
    viscircles([0,0],length(FRETdepth)+1,'color',[1,1,1],'LineWidth',35);
    hold off;
    axis equal
    axis([0 length(FRETdepth)+1 0 length(FRETdepth)+1]);
    title(['FRETretractionLineNormMean' ' t=' num2str(i17)],'FontSize',20);
    annotation('textbox',[.55,.8,.1,.1],'string',['FRET retraction t=' num2str(i17)],'FontSize',15,'FontWeight','bold');
    windowMax = gcf; windowMax.WindowState = 'Maximize';
    pause(.1);
    M(i17+halfT) = getframe();
    close all
end
v = VideoWriter([savepath 'FRETretractionLineNormMeanNorm.avi']);
v.FrameRate = 2;
open(v)
writeVideo(v,M)
close(v)

for i17 = timerange-n
    figure(1)
    hold on;
    for i18 = FRETdepth
        if ~isnan(FRETretractionMaxVLineNormMeanNorm(i18,i17+halfT))
            viscircles([0,0],length(FRETdepth)-i18+1,'color',[0,FRETretractionMaxVLineNormMeanNorm(i18,i17+halfT),0],'LineWidth',30);
        end
    end
    viscircles([0,0],length(FRETdepth)+1,'color',[1,1,1],'LineWidth',35);
    hold off;
    axis equal
    axis([0 length(FRETdepth)+1 0 length(FRETdepth)+1]);
    title(['FRETretractionMaxVLineNormMean' ' t=' num2str(i17)],'FontSize',20);
    annotation('textbox',[.55,.8,.1,.1],'string',['FRET retractionMaxV t=' num2str(i17)],'FontSize',15,'FontWeight','bold');
    windowMax = gcf; windowMax.WindowState = 'Maximize';
    pause(.1);
    M(i17+halfT) = getframe();
    close all
end
v = VideoWriter([savepath 'FRETretractionMaxVLineNormMeanNorm.avi']);
v.FrameRate = 2;
open(v)
writeVideo(v,M)
close(v)



figure(2)
% imagesc(FORCEprotrusionMean,[50,1000]);colormap jet;title(colorbar,'Pa');
imagesc(FORCEprotrusionMean);colormap jet;title(colorbar,'Pa');
xticks(1:6:length(timerange));xticklabels({'-120','-60','0','60','120'});yticks(1:5:length(FRETdepth));yticklabels({'1','6','11'});xlabel('t(s)');ylabel('depth(\mum)');
set(gca,'FontSize',16,'FontWeight','bold');
% title('FORCEprotrusionMean','FontSize',30);
windowMax = gcf; windowMax.WindowState = 'Maximize';
export_fig(gcf,[savepath 'FORCEprotrusionMean_' 'Distance' num2str(a) 'Valuerange' num2str(b) 'FRETdepth1to15FORCEdepth1to15.fig']);
export_fig(gcf,[savepath 'FORCEprotrusionMean_' 'Distance' num2str(a) 'Valuerange' num2str(b) 'FRETdepth1to15FORCEdepth1to15.png']);
close(2)

figure(2)
% imagesc(FORCEprotrusionMaxVMean,[50,1000]);colormap jet;title(colorbar,'Pa');
imagesc(FORCEprotrusionMaxVMean);colormap jet;title(colorbar,'Pa');
xticks(1:6:length(timerange));xticklabels({'-120','-60','0','60','120'});yticks(1:5:length(FRETdepth));yticklabels({'1','6','11'});xlabel('t(s)');ylabel('depth(\mum)');
set(gca,'FontSize',16,'FontWeight','bold');
% title('FORCEprotrusionMaxVMean','FontSize',30);
windowMax = gcf; windowMax.WindowState = 'Maximize';
export_fig(gcf,[savepath 'FORCEprotrusionMaxVMean_' 'Distance' num2str(a) 'Valuerange' num2str(b) 'FRETdepth1to15FORCEdepth1to15.fig']);
export_fig(gcf,[savepath 'FORCEprotrusionMaxVMean_' 'Distance' num2str(a) 'Valuerange' num2str(b) 'FRETdepth1to15FORCEdepth1to15.png']);
close(2)

figure(2)
% imagesc(FRETprotrusionMean,[1300,1900]);colormap jet;colorbar;
imagesc(FRETprotrusionMean);colormap jet;colorbar;
xticks(1:6:length(timerange));xticklabels({'-120','-60','0','60','120'});yticks(1:5:length(FRETdepth));yticklabels({'1','6','11'});xlabel('t(s)');ylabel('depth(\mum)');
set(gca,'FontSize',16,'FontWeight','bold');
% title('FRETprotrusionMean','FontSize',30);
windowMax = gcf; windowMax.WindowState = 'Maximize';
export_fig(gcf,[savepath 'FRETprotrusionMean_' 'Distance' num2str(a) 'Valuerange' num2str(b) 'FRETdepth1to15FORCEdepth1to15.fig']);
export_fig(gcf,[savepath 'FRETprotrusionMean_' 'Distance' num2str(a) 'Valuerange' num2str(b) 'FRETdepth1to15FORCEdepth1to15.png']);
close(2)

figure(2)
% imagesc(FRETprotrusionMaxVMean,[1500,1900]);colormap jet;colorbar;
imagesc(FRETprotrusionMaxVMean);colormap jet;colorbar;
xticks(1:6:length(timerange));xticklabels({'-120','-60','0','60','120'});yticks(1:5:length(FRETdepth));yticklabels({'1','6','11'});xlabel('t(s)');ylabel('depth(\mum)');
set(gca,'FontSize',16,'FontWeight','bold');
% title('FRETprotrusionMaxVMean','FontSize',30);
windowMax = gcf; windowMax.WindowState = 'Maximize';
export_fig(gcf,[savepath 'FRETprotrusionMaxVMean_' 'Distance' num2str(a) 'Valuerange' num2str(b) 'FRETdepth1to15FORCEdepth1to15.fig']);
export_fig(gcf,[savepath 'FRETprotrusionMaxVMean_' 'Distance' num2str(a) 'Valuerange' num2str(b) 'FRETdepth1to15FORCEdepth1to15.png']);
close(2)

figure(2)
% imagesc(FORCEretractionMean,[200,900]);colormap jet;title(colorbar,'Pa');
imagesc(FORCEretractionMean);colormap jet;title(colorbar,'Pa');
xticks(1:6:length(timerange));xticklabels({'-120','-60','0','60','120'});yticks(1:5:length(FRETdepth));yticklabels({'1','6','11'});xlabel('t(s)');ylabel('depth(\mum)');
set(gca,'FontSize',16,'FontWeight','bold');
% title('FORCEretractionMean','FontSize',30);
windowMax = gcf; windowMax.WindowState = 'Maximize';
export_fig(gcf,[savepath 'FORCEretractionMean_' 'Distance' num2str(a) 'Valuerange' num2str(b) 'FRETdepth1to15FORCEdepth1to15.fig']);
export_fig(gcf,[savepath 'FORCEretractionMean_' 'Distance' num2str(a) 'Valuerange' num2str(b) 'FRETdepth1to15FORCEdepth1to15.png']);
close(2)

figure(2)
% imagesc(FORCEretractionMaxVMean,[200,900]);colormap jet;title(colorbar,'Pa');
imagesc(FORCEretractionMaxVMean);colormap jet;title(colorbar,'Pa');
xticks(1:6:length(timerange));xticklabels({'-120','-60','0','60','120'});yticks(1:5:length(FRETdepth));yticklabels({'1','6','11'});xlabel('t(s)');ylabel('depth(\mum)');
set(gca,'FontSize',16,'FontWeight','bold');
% title('FORCEretractionMaxVMean','FontSize',30);
windowMax = gcf; windowMax.WindowState = 'Maximize';
export_fig(gcf,[savepath 'FORCEretractionMaxVMean_' 'Distance' num2str(a) 'Valuerange' num2str(b) 'FRETdepth1to15FORCEdepth1to15.fig']);
export_fig(gcf,[savepath 'FORCEretractionMaxVMean_' 'Distance' num2str(a) 'Valuerange' num2str(b) 'FRETdepth1to15FORCEdepth1to15.png']);
close(2)

figure(2)
% imagesc(FRETretractionMean,[1400,2000]);colormap jet;colorbar;
imagesc(FRETretractionMean);colormap jet;colorbar;
xticks(1:6:length(timerange));xticklabels({'-120','-60','0','60','120'});yticks(1:5:length(FRETdepth));yticklabels({'1','6','11'});xlabel('t(s)');ylabel('depth(\mum)');
set(gca,'FontSize',16,'FontWeight','bold');
% title('FRETretractionMean','FontSize',30);
windowMax = gcf; windowMax.WindowState = 'Maximize';
export_fig(gcf,[savepath 'FRETretractionMean_' 'Distance' num2str(a) 'Valuerange' num2str(b) 'FRETdepth1to15FORCEdepth1to15.fig']);
export_fig(gcf,[savepath 'FRETretractionMean_' 'Distance' num2str(a) 'Valuerange' num2str(b) 'FRETdepth1to15FORCEdepth1to15.png']);
close(2)

figure(2)
% imagesc(FRETretractionMaxVMean,[1400,1900]);colormap jet;colorbar;
imagesc(FRETretractionMaxVMean);colormap jet;colorbar;
xticks(1:6:length(timerange));xticklabels({'-120','-60','0','60','120'});yticks(1:5:length(FRETdepth));yticklabels({'1','6','11'});xlabel('t(s)');ylabel('depth(\mum)');
set(gca,'FontSize',16,'FontWeight','bold');
% title('FRETretractionMaxVMean','FontSize',30);
windowMax = gcf; windowMax.WindowState = 'Maximize';
export_fig(gcf,[savepath 'FRETretractionMaxVMean_' 'Distance' num2str(a) 'Valuerange' num2str(b) 'FRETdepth1to15FORCEdepth1to15.fig']);
export_fig(gcf,[savepath 'FRETretractionMaxVMean_' 'Distance' num2str(a) 'Valuerange' num2str(b) 'FRETdepth1to15FORCEdepth1to15.png']);
close(2)




figure(2)
% imagesc(FORCEprotrusionLineNormMean,[-0.3,0.3]);colormap jet;colorbar;
imagesc(FORCEprotrusionLineNormMean);colormap jet;colorbar;
xticks(1:6:length(timerange));xticklabels({'-120','-60','0','60','120'});yticks(1:5:length(FRETdepth));yticklabels({'1','6','11'});xlabel('t(s)');ylabel('depth(\mum)');
set(gca,'FontSize',16,'FontWeight','bold');
% title('FORCEprotrusionLineNormMean','FontSize',30);
windowMax = gcf; windowMax.WindowState = 'Maximize';
export_fig(gcf,[savepath 'FORCEprotrusionLineNormMean_' 'Distance' num2str(a) 'Valuerange' num2str(b) 'FRETdepth1to15FORCEdepth1to15.fig']);
export_fig(gcf,[savepath 'FORCEprotrusionLineNormMean_' 'Distance' num2str(a) 'Valuerange' num2str(b) 'FRETdepth1to15FORCEdepth1to15.png']);
close(2)

figure(2)
% imagesc(FORCEprotrusionMaxVLineNormMean,[-0.3,0.3]);colormap jet;colorbar;
imagesc(FORCEprotrusionMaxVLineNormMean);colormap jet;colorbar;
xticks(1:6:length(timerange));xticklabels({'-120','-60','0','60','120'});yticks(1:5:length(FRETdepth));yticklabels({'1','6','11'});xlabel('t(s)');ylabel('depth(\mum)');
set(gca,'FontSize',16,'FontWeight','bold');
% title('FORCEprotrusionMaxVLineNormMean','FontSize',30);
windowMax = gcf; windowMax.WindowState = 'Maximize';
export_fig(gcf,[savepath 'FORCEprotrusionMaxVLineNormMean_' 'Distance' num2str(a) 'Valuerange' num2str(b) 'FRETdepth1to15FORCEdepth1to15.fig']);
export_fig(gcf,[savepath 'FORCEprotrusionMaxVLineNormMean_' 'Distance' num2str(a) 'Valuerange' num2str(b) 'FRETdepth1to15FORCEdepth1to15.png']);
close(2)

figure(2)
% imagesc(FRETprotrusionLineNormMean,[-0.5,0.1]);colormap jet;colorbar;
imagesc(FRETprotrusionLineNormMean);colormap jet;colorbar;
xticks(1:6:length(timerange));xticklabels({'-120','-60','0','60','120'});yticks(1:5:length(FRETdepth));yticklabels({'1','6','11'});xlabel('t(s)');ylabel('depth(\mum)');
set(gca,'FontSize',16,'FontWeight','bold');
% title('FRETprotrusionLineNormMean','FontSize',30);
windowMax = gcf; windowMax.WindowState = 'Maximize';
export_fig(gcf,[savepath 'FRETprotrusionLineNormMean_' 'Distance' num2str(a) 'Valuerange' num2str(b) 'FRETdepth1to15FORCEdepth1to15.fig']);
export_fig(gcf,[savepath 'FRETprotrusionLineNormMean_' 'Distance' num2str(a) 'Valuerange' num2str(b) 'FRETdepth1to15FORCEdepth1to15.png']);
close(2)

figure(2)
% imagesc(FRETprotrusionMaxVLineNormMean,[-0.3,0.3]);colormap jet;colorbar;
imagesc(FRETprotrusionMaxVLineNormMean);colormap jet;colorbar;
xticks(1:6:length(timerange));xticklabels({'-120','-60','0','60','120'});yticks(1:5:length(FRETdepth));yticklabels({'1','6','11'});xlabel('t(s)');ylabel('depth(\mum)');
set(gca,'FontSize',16,'FontWeight','bold');
% title('FRETprotrusionMaxVLineNormMean','FontSize',30);
windowMax = gcf; windowMax.WindowState = 'Maximize';
export_fig(gcf,[savepath 'FRETprotrusionMaxVLineNormMean_' 'Distance' num2str(a) 'Valuerange' num2str(b) 'FRETdepth1to15FORCEdepth1to15.fig']);
export_fig(gcf,[savepath 'FRETprotrusionMaxVLineNormMean_' 'Distance' num2str(a) 'Valuerange' num2str(b) 'FRETdepth1to15FORCEdepth1to15.png']);
close(2)

figure(2)
% imagesc(FORCEretractionLineNormMean,[-0.3,0.3]);colormap jet;colorbar;
imagesc(FORCEretractionLineNormMean);colormap jet;colorbar;
xticks(1:6:length(timerange));xticklabels({'-120','-60','0','60','120'});yticks(1:5:length(FRETdepth));yticklabels({'1','6','11'});xlabel('t(s)');ylabel('depth(\mum)');
set(gca,'FontSize',16,'FontWeight','bold');
% title('FORCEretractionLineNormMean','FontSize',30);
windowMax = gcf; windowMax.WindowState = 'Maximize';
export_fig(gcf,[savepath 'FORCEretractionLineNormMean_' 'Distance' num2str(a) 'Valuerange' num2str(b) 'FRETdepth1to15FORCEdepth1to15.fig']);
export_fig(gcf,[savepath 'FORCEretractionLineNormMean_' 'Distance' num2str(a) 'Valuerange' num2str(b) 'FRETdepth1to15FORCEdepth1to15.png']);
close(2)

figure(2)
% imagesc(FORCEretractionMaxVLineNormMean,[-0.2,0.2]);colormap jet;colorbar;
imagesc(FORCEretractionMaxVLineNormMean);colormap jet;colorbar;
xticks(1:6:length(timerange));xticklabels({'-120','-60','0','60','120'});yticks(1:5:length(FRETdepth));yticklabels({'1','6','11'});xlabel('t(s)');ylabel('depth(\mum)');
set(gca,'FontSize',16,'FontWeight','bold');
% title('FORCEretractionMaxVLineNormMean','FontSize',30);
windowMax = gcf; windowMax.WindowState = 'Maximize';
export_fig(gcf,[savepath 'FORCEretractionMaxVLineNormMean_' 'Distance' num2str(a) 'Valuerange' num2str(b) 'FRETdepth1to15FORCEdepth1to15.fig']);
export_fig(gcf,[savepath 'FORCEretractionMaxVLineNormMean_' 'Distance' num2str(a) 'Valuerange' num2str(b) 'FRETdepth1to15FORCEdepth1to15.png']);
close(2)

figure(2)
% imagesc(FRETretractionLineNormMean,[-0.1,0.7]);colormap jet;colorbar;
imagesc(FRETretractionLineNormMean);colormap jet;colorbar;
xticks(1:6:length(timerange));xticklabels({'-120','-60','0','60','120'});yticks(1:5:length(FRETdepth));yticklabels({'1','6','11'});xlabel('t(s)');ylabel('depth(\mum)');
set(gca,'FontSize',16,'FontWeight','bold');
% title('FRETretractionLineNormMean','FontSize',30);
windowMax = gcf; windowMax.WindowState = 'Maximize';
export_fig(gcf,[savepath 'FRETretractionLineNormMean_' 'Distance' num2str(a) 'Valuerange' num2str(b) 'FRETdepth1to15FORCEdepth1to15.fig']);
export_fig(gcf,[savepath 'FRETretractionLineNormMean_' 'Distance' num2str(a) 'Valuerange' num2str(b) 'FRETdepth1to15FORCEdepth1to15.png']);
close(2)

figure(2)
% imagesc(FRETretractionMaxVLineNormMean,[-0.25,0.25]);colormap jet;colorbar;
imagesc(FRETretractionMaxVLineNormMean);colormap jet;colorbar;
xticks(1:6:length(timerange));xticklabels({'-120','-60','0','60','120'});yticks(1:5:length(FRETdepth));yticklabels({'1','6','11'});xlabel('t(s)');ylabel('depth(\mum)');
set(gca,'FontSize',16,'FontWeight','bold');
% title('FRETretractionMaxVLineNormMean','FontSize',30);
windowMax = gcf; windowMax.WindowState = 'Maximize';
export_fig(gcf,[savepath 'FRETretractionMaxVLineNormMean_' 'Distance' num2str(a) 'Valuerange' num2str(b) 'FRETdepth1to15FORCEdepth1to15.fig']);
export_fig(gcf,[savepath 'FRETretractionMaxVLineNormMean_' 'Distance' num2str(a) 'Valuerange' num2str(b) 'FRETdepth1to15FORCEdepth1to15.png']);
close(2)




close all

save([savepath 'Distance' num2str(a) 'Valuerange' num2str(b) 'FRETdepth1to15FORCEdepth1to15final.mat']);



['Distance' num2str(a) 'Valuerange' num2str(b) 'FRETdepth1to15FORCEdepth1to15__runtime(min):' num2str(toc/60)]

