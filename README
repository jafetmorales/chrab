Welcome to Chrab

Chrab is a Python3 program that measures chromatic aberration.

How to use:

**Step 1:
Generate a chromatically aberrated image using Kromo by executing the following command:
python3 kromo.py --verbose --strength 1 --noblur --out flower_chrabed.jpg ./samples/flower.jpg

**Step 2:
Measure chromatic aberration for your image as follows:
python3 chrab.py -i ./flower_chrabed.jpg -c red
The output will show the coordinates of the center for the aberration and the alpha parameter
for the aberration.


**How to get more accurate results:
To get more accurate results you can adjust the side, delta, and K variables in chrab.py:

-SIDE: The side variable in chrab.py determines the dimensions of the sampling window that is used to measure mutual information. Ideally,
you would have a high number for side but this is very computationally expensive because the number of pixel positions to be warped is side^2.
Program could generate an error if the sampling window goes outside of the image boundaries.

-DELTA: The delta variable in chrab.py determines the stride when moving the sampling window that is used to measure mutual information. Ideally,
you would have a small delta so that you scan more places inside the image, however, this requires much more computations.
Program could generate an error if the sampling window goes outside of the image boundaries.

-K: The K variable in chrab.py determines in how many places you will measure mutual information. Ideally you would sample the entire image,
but this is very computationally expensive because the number of samples locations will be 4K^2.
Program could generate an error if the sampling window goes outside of the image boundaries.
