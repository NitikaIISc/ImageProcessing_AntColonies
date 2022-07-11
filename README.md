# ImageProcessing_AntColonies
Protocol for data collection and image analysis to understand the interplay of nest architecture and colony size on collective food distribution in Argentine ant colonies. Code written by Alon Oyler-Yaniv and modified by Nitika Sharma and Julie Miller
The location of ants themselves can also be extracted at different timepoints within the colony using the red channel in the images as the pixels below mean - x standard deviation as pixels with ants (Image 3). 
 
See AntsExperimentProtocol.pdf
Image 3: The location of ants extracted using the red channel of images can be viewed in purple (left) and the fluorescent fool in the abdomens of ants as well as being exchanged via trophallaxis can be viewed in white pixels (right).
NOTE: I identified some issues that can be troubleshot by playing around with the standard deviation in the threshold being set in the MATLAB code to identify ants and food. Troubleshooting is a trade-off between false positives and false negatives and varies across videos (see examples in Image 4 below). 
 
Image 4: Some variation in video quality and accompanying issues that needed troubleshooting.
We dealt with some of the issues like lighting reflection creating blind spots in the images by incorporating a code that will allow hand-drawn polygons with blind spots to be excluded from the image analysis, so they are not mistaken for fluorescent food. The code also allows to differentiate between different nest shaped that were filmed by asking for the number of wells in a selected video, which would be 1 for circular nests, 4 for clover nests and 9 for the big 9-chambered nests and can be hand-drawn in MATLAB code (Ver5_20210227_NitikaTroubleshotProcessTifFolder.m) when prompted.
To avoid slow processing and huge files with coordinates for negative space, I extracted only the pixels that contained ants and food in the progressive images which grow exponentially as the food gets distributed within the colony. 
