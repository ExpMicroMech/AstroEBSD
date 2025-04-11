function tiledPattern = TileQuadPatterns(quadrantData)

quadrantSize = 256;
quadrantpatternSize = quadrantSize ^ 2;
fullRawPatternSize = quadrantSize ^ 2 * 4;

edgeoffset = 35;

rawPatternSize = 549;
patternSize = 479;
magicNumber = patternSize * patternSize;

max_quadrantData = max(quadrantData(:));
 
ch0 = reshape(quadrantData(1:quadrantpatternSize), quadrantSize, quadrantSize)';
ch1 = reshape(quadrantData(quadrantpatternSize+1:2*quadrantpatternSize), quadrantSize, quadrantSize)';
ch2 = reshape(quadrantData(2*quadrantpatternSize+1:3*quadrantpatternSize), quadrantSize, quadrantSize)';
ch3 = reshape(quadrantData(3*quadrantpatternSize+1:4*quadrantpatternSize), quadrantSize, quadrantSize)';

tiledPattern = zeros(rawPatternSize) + max_quadrantData;

tiledPattern(1:256, 36:256+35) = rot90(ch3, 1);
tiledPattern(36:256+35, 256+2+36:256+2+35+256) = rot90(ch0, 3);
tiledPattern(256+2+1:512+2, 1:256) = rot90(ch2, 2);
tiledPattern(256+2+36:512+2+35, 256+2+1:512+2) = rot90(ch1, 3);

tiledPattern = tiledPattern(edgeoffset+1:rawPatternSize-edgeoffset, edgeoffset+1:rawPatternSize-edgeoffset);
end