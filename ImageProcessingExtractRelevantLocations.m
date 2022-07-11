%% Load all images in the folder, convert to RGB

%baseDir = 'F:/Nitika_architecture/20200627_col7_clover_2306/post-food/_1/'
baseDir = 'Z:/Nitika/NewlyTransfer/Ant3_col7to14/Julie and Nitika_Col7to14/Temescal_7_small_circle_09102020_postfluor_1'
%baseDir = 'F:/Nitika_architecture/20200729_col11_circle_2306/post-food/' %parent directory where all the movies are
%baseDir = 'C:/Users/nitik/Box Sync/Manuscript6_Julie/JuliesAnts-master'
myFolder = uigetdir(baseDir);

if ~isdir(myFolder)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
    uiwait(warndlg(errorMessage));
    return;
end
%find all *.tif files
filePattern = fullfile(myFolder, '*.tif');
tifFiles = dir(filePattern);

%make placeholder
% stack = cell(numel(tifFiles),1);
% 
% %load and convert all images
% for k = 1:length(tifFiles)
%     baseFileName = tifFiles(k).name;
%     fullFileName = fullfile(myFolder, baseFileName);
%     fprintf(1, 'Now reading %s\n', fullFileName);
%     imageArray = mosaicToRGB(single(imread(fullFileName))/2^16);
%     stack{k} = imageArray;
% end
% stack = cat(4,stack{:});

%% Look at the first timeframe. Create masks for the different wells

k=17000;
%clover = 1533;

fullFileName = fullfile(myFolder, tifFiles(k).name);
fprintf(1, 'Now reading %s\n', fullFileName);
im = mosaicToRGB(single(imread(fullFileName))/2^16);



hFig = figure(332);
hFig.NumberTitle = 'off';
hFig.Name = 'How many wells'
imagesc(im); shg
prompt = 'How many wells? ';
nWells = str2double(inputdlg(prompt));

%BW will be a mask (1 inside the wells, 0 outside)
BW = zeros(size(im,1), size(im,2));

for i=1:nWells
    hFig.Name = ['Mask well ' num2str(i) ' of ' num2str(nWells)]
    if i==1
        [bw1, el] = roicirclecrop(im);
        
    else
        [bw1, el] = roicirclecrop(bsxfun(@times, im,~BW),'Position', el.getPosition);
    end
    BW = BW + bw1;
    BW = BW>0;
end
imagesc(bsxfun(@times, im,~BW)) %viewing manually cropped well to consider for image analysis


%% Exclude regions

hFig = figure(332);
hFig.NumberTitle = 'off';
hFig.Name = 'Find regions to exclude'
imagesc(bsxfun(@times, im,BW)); shg
prompt = 'Are there regions to exclude?'

regionstoexclude = menu(prompt,'Yes','No');

%As long as there is at least one region to exclude because of reflection
%of light on plexiglass, keep selecting
while regionstoexclude==1  
    [bw1, el] = roipolycrop(bsxfun(@times, im,BW));
        BW = BW - bw1;
    BW = BW>0;
    imagesc(bsxfun(@times, im,BW)); shg
    regionstoexclude = menu(prompt,'Yes','No');
end

hFig.Name = 'Final masked image'

imagesc(bsxfun(@times, im,BW))


%% Find threshold for what is considered "food"
%im is now an RGB image, the green channel is saturated, the blue channel
%has reflections from the light source. The red channels is ok. Ideally
%there should be some hardware fixes for this in the future.

%reloading. You may want to sometimes use another timepoint as your ref so
%i'm leaving it here
k=23999
% spiral = 16108;

%%Now reading images inside this folder
fullFileName = fullfile(myFolder, tifFiles(k).name);
fprintf(1, 'Now reading %s\n', fullFileName);
im = mosaicToRGB(single(imread(fullFileName))/2^16);

I = im(:,:,1); %red channel

hFig = figure(332);
hFig.NumberTitle = 'off';
hFig.Name = 'Select region inside the wells with NO ants' %%setting negative space in image

imagesc(I)
axialROI=impoly;
axialROImask=createMask(axialROI);

%taking anything above mean+x std in the red channel as food
threshFood = mean(I(axialROImask))+0.25*std(I(axialROImask));
imshowpair(im(:,:,1), I.*(I>threshFood));
hFig.Name = 'White is food'

%taking anything below mean-x std in the red channel as Ants
hFig2 = figure(333);
hFig2.NumberTitle = 'off';

threshAmIAnt = mean(I(axialROImask))-4.15*std(I(axialROImask));
imshowpair(I, ((I.*BW)<threshAmIAnt-(~BW)));
hFig2.Name = 'Purple is Ants'




%% Now process whole stack

totalFood = nan(1,numel(tifFiles));
totalAnts = nan(1,numel(tifFiles));

fullFileName = fullfile(myFolder, tifFiles(1).name);
fprintf(1, 'Now reading %s\n', fullFileName);
im = mosaicToRGB(single(imread(fullFileName))/2^16);
imsize = size(im);

ImgBin = 500
FoodLocs = [];
AntsLocs = [];
%for k=1:10;
for k=1:numel(tifFiles);
    try
        
        fullFileName = fullfile(myFolder, tifFiles(k).name);
        fprintf(1, 'Now reading %s\n', fullFileName);
        im = mosaicToRGB(single(imread(fullFileName))/2^16);
        %bsxfun(function, A,B); @times is .* element by element
        %multiplication
        FoodMat=bsxfun(@times,(im(:,:,1)-median(I(axialROImask))).*(im(:,:,1)>threshFood), BW); %.*multiplying two matrices?
        AntsMat=bsxfun(@minus, bsxfun(@times, squeeze(im(:,:,1)), BW)<threshAmIAnt,~BW);
        totalFood(k) = squeeze(sum(sum(FoodMat(k))));
        totalAnts(k) = squeeze(sum(sum(AntsMat(k))));
        %add k as a column
        [rows, columns, values] = find(FoodMat(:,:));
        %FoodLocs{k} = [FoodLocs; rows, columns, values]
        tempk = k*ones(length(rows),1);
        FoodLocs = [FoodLocs;rows, columns, values,tempk];
        %%if b = mod (k,500) = if b = 0
        %% writematrix(paste(k+FoodLocs))<-for every b = 0
        [rowsA, columnsA, valuesA] = find(AntsMat(:,:));
        %FoodLocs{k} = [FoodLocs; rows, columns, values]
        tempkA = k*ones(length(rowsA),1);
        AntsLocs = [AntsLocs;rowsA, columnsA, valuesA,tempkA];
        %[row,col,v] = find(___) also returns vector v, which contains the nonzero elements of X.
        %XYFood(:,:,k) = [rows, columns, values]
        
        %%%Nitika edit:
        %for every 5 images' food and ant matrices appended, export it to
        %the folder by naming the matrix appropriately
        
        if mod(k, ImgBin) == 0
            
          writematrix(FoodLocs,['F:/Julie and Nitika Videos/Kerckhoff_15_big_circle_11082020_postfluor_1/Processed_locationMatrices/' 'FoodLocs' num2str(k) '.csv'])
          writematrix(AntsLocs,['F:/Julie and Nitika Videos/Kerckhoff_15_big_circle_11082020_postfluor_1/Processed_locationMatrices/' 'AntsLocs' num2str(k) '.csv'])
        %empty after every 500 images so that memory doesn't exhaust
          FoodLocs = [];
          AntsLocs = [];
        end
        
    catch
        warning(['couldn`t load file ' fullFileName] )
    end
end

totalFood(isnan(totalFood))=[];
totalAnts(isnan(totalAnts))=[];


%% Plot results
set(0,'DefaultTextInterpreter', 'tex')
set(0, 'DefaultAxesFontName', 'Arial')
set(0, 'DefaultAxesFontSize', 20)
set(0, 'DefaultUIControlFontName', 'Arial')
set(0,'defaulttextfontname','Arial');
set(0,'defaulttextfontsize',22);
set(groot,'defaultFigureColor','w')
set(groot,'defaultAxesColor','w')
set(groot,'DefaultLineMarkerSize',3)
set(groot,'defaultAxesTickLength',[0.03 0.01])
set(groot,'defaultLineLineWidth',2)

yyaxis left
plot(1:length(totalFood), totalFood)
ylabel('total Food(a.u.)')
yyaxis right
plot(1:length(totalFood), totalAnts); shg
xlabel('time(frames)')
ylabel('total Ants(a.u.)')


%viewing food and ants
for i = 17990:17999 imshowpair(AntsMat(:,:), FoodMat(:,:)); shg; pause(1); end

className = class(FoodMat);
disp(className)
size(FoodMat)


%%Exporting for use in R


%save as csv
writematrix(FoodLocs,'Z:/Nitika/NewlyTransfer/Ant3_col7to14/Julie and Nitika_Col7to14/Temescal_7_small_circle_09102020_postfluor_1/Processed/FoodLocsAll_excluded.csv')
writematrix(AntsLocs,'Z:/Nitika/NewlyTransfer/Ant3_col7to14/Julie and Nitika_Col7to14/Temescal_7_small_circle_09102020_postfluor_1/Processed/AntsLocsAll_excluded.csv')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% Write h5file
h5create('myfile.h5','/DS1',[10 20])
mydata = rand(10,20);
h5write('myfile.h5', '/DS1', mydata)

h5create('FoodH5.h5','/FoodMatrix',size(FoodMat))
h5write('FoodH5.h5','/FoodMatrix',FoodMat)
h5disp('FoodH5.h5')


h5create('FoodLocs.h5','/FoodLocs',size(FoodLocs,2))
h5write('FoodLocs.h5','/FoodLocs',FoodLocs{1:size(FoodLocs,2)})
h5disp('FoodLocs.h5')
%Vertically combine cell arrays %check out:https://www.mathworks.com/matlabcentral/answers/253988-concatenate-unequal-sized-arrays-produced-in-a-loop
AllFoodLocs = FoodLocs{1:size(FoodLocs,2)}


%save('Food500.mat', 'FoodMat', '-v7.3')
save('Food500V6.mat', 'FoodMat', '-v6')

%FoodMat and AntMat are called "Multidimensional arrays"

%save as matlab file
save Food500.mat FoodMat


% Credit to cyclist for this
C = permute(A,[1 3 2]);
C = reshape(C,[],size(A,2),1)

size(FoodMat)
TwoDFoodMat = reshape(FoodMat,[size(FoodMat,1)*size(FoodMat,3),size(FoodMat, 2)])
save('TwoDFoodMat.mat', 'TwoDFoodMat', '-v6')
writematrix(TwoDFoodMat,'TwoDFoodMat.csv') 


a = FoodMat

con = fopen('a.bin', 'w');
fwrite(con, a * 0.01, 'float64')
fclose(con)

a * 0.01





A = (reshape(FoodMat,768,1024,[]))
size(A)
%clear;clc;
M= FoodMat;
%[1000000,12,2,2];
dlmwrite('a.txt',M); % save M to file--a.txt
type a.txt; % print content in a.txt
M = dlmread('a.txt'); % load content of a.txt to M and then you will have 'M=[1000000,12,2,2]'



% append csv files as each matrix
for i = 1:(numel(tifFiles)-1);
% Write the table to a CSV file
writematrix(FoodMat(:,:,i),'a.csv');
writematrix(FoodMat(:,:,i+1),'a.csv','WriteMode','append'); end