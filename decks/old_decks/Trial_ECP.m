clear
close all
home

imagefolder='C:\Users\Public\Desktop\SharedData\Haroon\2025\0806_SAECP_Alignments';
image1='20keV_3.2nA_atCircles';
image2='20keV_3.2nA_Ref';

image_1_full=fullfile(imagefolder,[image1 '.tif']);
image_2_full=fullfile(imagefolder,[image2 '.tif']);

i1=imread(image_1_full);
i2=imread(image_2_full);

info1 = imfinfo(image_1_full);
info2 = imfinfo(image_2_full);

y_height=2/3*size(i1,2);

i1_image_only=i1(1:y_height,:);

figure; imagesc(i1_image_only); axis image; axis tight; colormap('gray');
info1 = imfinfo(image_1_full);
i1_text=info1.UnknownTags.Value;
i1_text_pairs=regexp(i1_text, '[\f\n\r]', 'split');
nonEmptyCells = ~cellfun('isempty', i1_text_pairs);
i1_text_pairs_full=i1_text_pairs(nonEmptyCells);

%find the text

% contains_text = cellfun(@(x)strfind(x,'StageRawR'), i1_text_pairs);
[text_retrieved] = fReadTFSHeaderPair('StageRawR',i1_text_pairs_full);

%%
figure;
subplot(1,3,1);
imagesc(i1)
title(image1)
axis equal; axis image; colormap('gray');

subplot(1,3,2);
imagesc(i2)
title(image2)
axis equal; axis image; colormap('gray');

% Display the difference between the two images
subplot(1,3,3);
imagesc(abs(i1 - i2));
title('Difference Image');
axis equal; axis image; colormap('jet');