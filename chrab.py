import math
from PIL import Image 
import numpy as np
import scipy as sp
import scipy.interpolate
import sys
import warnings
import string
import sys, getopt

# LOAD THE JPEG IMAGE
def ldimage(inputFilePath):
    img = Image.open(inputFilePath) 
    return img
    
# CALCULATE MUTUAL INFORMATION
def mut_info(hgram):
    """ Mutual information for joint histogram
    """
    # Convert bins counts to probability values
    pxy = hgram / float(np.sum(hgram))
    px = np.sum(pxy, axis=1) # marginal for x over y
    py = np.sum(pxy, axis=0) # marginal for y over x
    px_py = px[:, None] * py[None, :] # Broadcast to multiply marginals
    # Now we can do the calculation using the pxy, px_py 2D arrays
    nzs = pxy > 0 # Only non-zero pxy values contribute to the sum
    return np.sum(pxy[nzs] * np.log(pxy[nzs] / px_py[nzs]))

# THIS CLASS CALCULATES AND STORES THE PARAMS FOR A SAMPLING WINDOW.
# IT MAKES IT EASIER TO TRANSFER THESE PARAMETERS BETWEEN FUNCTIONS.
class RegionParams:
    cx=0
  
def getRegionParams(img, cx,cy,side):
    top_row=cy-int(round(side/2))
    bottom_row=cy+int(round(side/2))
    left_col=cx-int(round(side/2))
    right_col=cx+int(round(side/2))
    regParams=RegionParams()
    regParams.cx=cx
    regParams.cy=cy
    regParams.top_row=top_row
    regParams.bottom_row=bottom_row
    regParams.left_col=left_col
    regParams.right_col=right_col
    return regParams
    
# THIS GETS A REGION FROM AN IMAGE USING THE PARAMETERS IN A RegionParams OBJECT
def getRegion(img,regionParams):
    region=img[regionParams.top_row:regionParams.bottom_row+1, regionParams.left_col:regionParams.right_col+1]
    return region

# THIS FUNCTION WARPS A REGION IN THE GREEEN CHANNEL. THE WARPING IS DONE BY MOVING PIXEL POSITIONS BY A GIVEN ALPHA VALUE
# AFTER WE MOVE THE PIXEL POSITIONS, THE PIXEL POSITIONS WILL FALL IN PLACES THAT ARE NOT INTEGERS AND WHICH DO NOT CORRESPONDS TO
# TO THE PIXEL POSITIONS OF A REGULAR IMAGE [0,0][0,1],...[WIDTH,HEIGHT]. THEREFORE, WE MUST USE INTERPOLATION TO INTERPOLATE THE PIXEL
# VALUES AT THE PIXEL POSITIONS OF A REGULAR IMAGE. THIS WILL ALLOW US TO CALCULATE MUTUAL INFORMATION LATER ON.
def warpRegion(region, regParams,alpha):
    x_nowarp, y_nowarp = np.mgrid[regParams.left_col:regParams.right_col+1, regParams.top_row:regParams.bottom_row+1]
    x_nowarp_flat=x_nowarp.flatten()
    y_nowarp_flat=y_nowarp.flatten()
    #HERE WE WARP THE PIXEL POSITIONS (BY WARPING THE INDEXES OF THE PIXELS)
    x_warped_flat = alpha*(x_nowarp_flat - regParams.cx)+regParams.cx;
    y_warped_flat = alpha*(y_nowarp_flat - regParams.cy)+regParams.cy;
    #VECTORIZE THE INDEXES SINCE THEY WERE IN A 2D ARRAY AND WE NEED THEM AS VECTORS
    flat_region=region.flatten()
    #MAKE NUMPY ABLE TO PRINT ARRAYS COMPLETELY
    np.set_printoptions(threshold=sys.maxsize)
    #METHOD 1 FOR INTERPOLATION (faster but less accurate because uses linear interpolation or first order spline)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        bSpline=scipy.interpolate.bisplrep(x_warped_flat,y_warped_flat,flat_region, w=None, xb=None, xe=None, yb=None, ye=None, kx=1, ky=1, task=0, s=None, eps=1e-16, tx=None, ty=None, full_output=0, nxest=None, nyest=None, quiet=1)
    flat_region_warped = [float(scipy.interpolate.bisplev(XX, YY, bSpline, dx=0, dy=0)) for XX,YY in zip(x_nowarp_flat,y_nowarp_flat)]
    #METHOD 2 FOR INTERPOLATION (takes too long but more accurate because uses cubic or third order spline)
    # interpol = sp.interpolate.interp2d(x_warped_flat,y_warped_flat,flat_region, 'cubic');
    # flat_region_warped = [float(interpol(XX,YY)) for XX,YY in zip(x_nowarp_flat,y_nowarp_flat)]
    #METHOD 3 FOR INTERPOLATION (don't do it. quintic interpolation will most likely take forever)
    
    # RESHAPE THE IMAGE TO BE A SQUARE AGAIN
    region_warped=np.reshape(flat_region_warped, region.shape) # C-like index ordering
    return region_warped
    
# THE MAIN
def main(argv):
    # TAKE IN THE PARAMS
    inputfile = ''
    channel = ''
    try:
      opts, args = getopt.getopt(argv,"hi:c:",["ifile=","channel="])
    except getopt.GetoptError:
      print('measure.py -i <inputfile> -c <red or blue>')
      sys.exit(2)
    for opt, arg in opts:
      if opt == '-h':
         print('measure.py -i <inputfile> -c <red or blue>')
         sys.exit()
      elif opt in ("-i", "--ifile"):
         inputfile = arg
      elif opt in ("-c", "--channel"):
         channel = arg
    print('Input file is ', inputfile)
    print('Channel is ', channel)
   
    np. set_printoptions(suppress=True)
    
    # LOAD IMAGE
    img=ldimage(inputfile)
    # SPLIT INTO BANDS
    r, g, b = img.split()
    
    #creating an array of pixels so that we can access them
    gg = np.array(g)

    # TAKE THE CHANNEL THE USER WANTS
    if(channel=='red'):
        cc=np.array(r)
    elif(channel=='blue'):
        cc=np.array(b)
    else:
        raise Exception("Please choose either red or blue for the channel for which you want to estimate chromatic aberration")
    
    # ADJUST THIS PARAMETERS AS INDICATED IN README IF YOU WANT MORE ACCURATE RESULTS. LARGE NUMBERS CAN MAKE IT COMPUTAIONALLY
    # EXPENSIVE DEPENDING ON THE COMPLEXITY OF THE PARAMETER.
    side=50
    delta=50
    K=6
    histBins=255
    
    # DETERMINE IMAGE CENTER, SO WE CAN MOVE THE SAMPLING WINDOW AROUND IT 
    imgCenterX=int(round(gg.shape[1]/2))
    imgCenterY=int(round(gg.shape[0]/2))

    aberrationEstimates=np.empty((0,4))
    print("Chromatic aberration estimates [cx,cy,alpha,mutual information]:")
    # MOVE THE SAMPLING WINDOW ACROSS THE IMAGE
    for cx in np.arange(imgCenterX-K*delta,imgCenterX+K*delta,delta):
        for cy in np.arange(imgCenterY-K*delta,imgCenterY+K*delta,delta):
            regParams=getRegionParams(gg, cx, cy, side)
            #GET THE REGIONS TO BE COMPARED (THEY WILL BE COMPARED AFTER ONE OF THEM IS WARPED LATER ON) FROM TWO CHANNELS 
            region_gg_unwarped=getRegion(gg,regParams)
            region_cc_unwarped=getRegion(cc,regParams)
            #CALCULATE MUTUAL INFORMATION AT DESIRED REGION BETWEEN UNWARPED RED CHANNEL AND WARPED GREEN CHANNEL FOR DIFFERENT WARPING SCALES (DIFFERENT ALPHAS) TO SEE HOW WARPED THEY ARE 
            for alpha in np.arange(.5, 1.5, 0.2):
                # DO THE WARPING AT NO WARP SPEED (IT'S COMPUTATIONALLY EXPENSIVE)
                region_gg_warped=warpRegion(region_gg_unwarped,regParams,alpha)
                #CALCULATE MUTUAL INFORMATION BETWEEN UNWARPED RED CHANNEL AND WARPED GREEN CHANNEL
                joinHist, x_edges, y_edges = np.histogram2d(region_cc_unwarped.ravel(),region_gg_warped.ravel(),bins=histBins)
                minfo=mut_info(joinHist)
                #APPEND TO THE ESTIMATES VECTOR 
                aberrationEstimates=np.append(aberrationEstimates,[[cx,cy,alpha,minfo]], axis=0)
                # PRINT TO AVOID DESPERATION OF USER
                print([cx,cy,alpha,minfo])

    # CHOOSE THE PARAMETERS THAT GIVE THE MAXIMUM MUTUAL INFORMATION
    maxMinfoIndex=np.argmax(aberrationEstimates[:,3],axis=0)
    aberrationBestEstimate=aberrationEstimates[maxMinfoIndex,:]
    
    print("The chromatic aberration estimate is:")
    print("Format: [cx, cy, alpha, mutual information]")
    print(aberrationBestEstimate)

    # FOR DEBUGGING PURPOSES WE CAN STORE THE WARPED AND UNWARPED SAMPLING WINDOWS. THE WARPED IMAGE SHOULD LOOK LIKE A ZOOM IN OF THE UNWARPED IMAGE, SINCE IT IS SCALED.
    imgRegion_warped = Image.fromarray(region_gg_warped).convert('RGB')
    imgRegion_warped.save("region_warped.jpeg")
    imgRegion=Image.fromarray(region_gg_unwarped).convert('RGB')
    imgRegion.save("region.jpeg")
   
if __name__ == "__main__":
   main(sys.argv[1:])



