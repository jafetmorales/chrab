import math
from PIL import Image 
import numpy as np
from decimal import Decimal
import scipy as sp
from scipy import interpolate
from scitools.std import ndgrid
from scipy import ogrid, sin, mgrid, ndimage, array


def ldimage():
    #load image
    global im
    # im = Image.open("/home/areej/Desktop/mandril_color.tif") 
    im = Image.open("./samples/flower_chromatic.jpg") 

def analyzeCA(mode, im):
    n_regions = 10;
    reg_size = [300, 300];
    overlap = 0.5;


    levels = 9;
    steps = 2;
    edge_width = 10;
    hist_sz = 128;

   # alpha_1 and alpha_2 are assumed to be between these values
    w_data = [0.9985, 1.0015];

    reg_list=[]

   #creating an array of pixels so that we can access them
    pix=im.load()

    #
    #Analyze full image
    if mode=='full':
        print("Doing a full analysis")
        # mx_shift is the third argument in 'full' mode
        mx_shift = n_regions;
            # [ydim,xdim,zdim]= size(im);
        ydim=im.size[0]
        xdim=im.size[1]
        zdim=3

        #print "Image dimensions: [ydim, xdim, zdim]= "+str(ydim)+' '+str(xdim)+' '+str(zdim)#str([ydim,xdim,zdim])
        print("Image dimensions: [ydim, xdim, zdim]= "+str([ydim,xdim,zdim]))


        global alpha_mx, alpha_mn
        alpha_mx = 1 + 4*mx_shift / math.sqrt( xdim*xdim + ydim*ydim );
        alpha_mn = 1.0/alpha_mx;

        print("alpha_mx= "+str(alpha_mx))
        print("alpha_mn= "+str(alpha_mn))
        #recompute alpha_1 and alpha_2 to be between
        #these new values
        w_data = [alpha_mn, alpha_mx];
        ew = edge_width;

        #take the image minus a ew-wide edge
        roi = [ew+1, xdim-ew, ew+1, ydim-ew];

        print("edge_width= "+str(ew))
        print("roi= "+str(roi))

        #Analyze blue to green chromatic aberration
        bg_params = parameterSearch( im, [3, 2], roi, ew, hist_sz, w_data);

        # Analyze red to green chromatic aberration
        rg_params = parameterSearch( im, [1, 2], roi, ew, hist_sz, w_data );
    elif mode=='reg':
      print("we should do a regional analysis here")

    else:
     print("unsupported call")

#def estimateCARegions( im, [3, 2], reg_list, settings ):
def parameterSearch( im, colour_space, roi, ew, hist_sz, w_data):

    #levels is number of iterations 
    levels = 8;
    steps = 2;

    #[ydim,xdim,zdim] = size(im);
    ydim=im.size[0]
    xdim=im.size[1]
    zdim= 3


    x_data = [1, xdim];
    y_data = [1, ydim];

    xlim = x_data;
    ylim = y_data;
    zlim = w_data;

#work out which of height and width is the bigger
    dim = max(xdim,ydim)

    print("The highest dimension is : "+str(dim))

#check that roi falls within expected boundries
    if ((roi[0] <= ew) or (roi[1] > xdim-ew) or (roi[2] <= ew) or (roi[3] > ydim-ew)):
        print("ROI is too close to image edges")
        return -1 # TODO: terminate here with an error
        #Get image regions

    source = im.split()
    Cfixed = source[2]
    Cwarp  = source[1]
    #[ydim,xdim,zdim] = size(im);
    ydimCwarp=Cwarp.size[0]
    xdimCwarp=Cwarp.size[1]
    print('xdimCwarp'+str(xdimCwarp))

    roi_pad = [roi[0]-ew, roi[1]+ew, roi[2]-ew, roi[3]+ew];
    for levels in range(1,8):
        #Guess at a center and then compute best warp
        #user defined function linear_space used to generate linearly spaced vectors
        x_coords = np.linspace(0,511,steps+2)
        y_coords = np.linspace(0,511,steps+2)
        z_coords = np.linspace(alpha_mn,alpha_mx,steps+2)
        step_x=(xlim[1]-xlim[0])/(steps+1)
        start_x=xlim[0]+step_x
        end_x=xlim[1]-step_x+0.5
        step_y=(ylim[1]-ylim[0])/(steps+1)
        start_y=ylim[0]+step_y
        end_y=ylim[1]-step_y+0.5
        step_z=(zlim[1]-zlim[0])/(steps+1)
        start_z=zlim[0]+step_z
        fudge_z=step_z/2.0
        end_z=zlim[1]-step_z+fudge_z
        #Do not include end points in search;
        centers_x, centers_y, warps= np.mgrid[start_x:end_x:step_x,start_y:end_y:step_y,start_z:end_z:step_z]
        centers_x=centers_x.flatten()
        centers_y=centers_y.flatten()
        warps=warps.flatten()
        mi = np.zeros(centers_x.size)

        for k in range(0,centers_x.size):
           cx = centers_x[k]
           cy = centers_y[k]
           wz = warps[k]       
           #Warp the region 
           temp_im = warpRegion(Cwarp, roi_pad, cx, cy, wz)
                #correlation
           mi[k] = np.corrcoef(Cfixed, temp_im)
               #Now pick the best quadrant
        v, max_ix = math.max(mi)
        ix, jx, kx = arrayInd(mi.size, max_ix);
        ##The coordinates of err are off by 1 from x_coords and y_coords because
        ##we did not include the end point
    xlim = x_coords([jx, jx+2]);
    ylim = y_coords([ix, ix+2]);
    zlim = z_coords([kx, kx+2]);

    cx = math.mean(xlim);
    cy = math.mean(ylim);
    wz = math.mean(zlim);

    print("x= "+str(cx))
    print("y= "+str(cy))
    print("z= "+str(wz))
    
def warpRegionOLDDDDD(Cwarp, roi_pad, cx, cy, wz):
#Unpack region indices
    sx, ex, sy, ey = roi_pad

    xramp, yramp = np.mgrid[sx:ex+1, sy:ey+1]

    xrampc = xramp - cx;
    yrampc = yramp - cy;
    xramp1 = 1/wz*xrampc;
    yramp1 = 1/wz*yrampc;
    xrampf = xrampc.flatten()
    yrampf = yrampc.flatten()
    xramp1f = xramp1.flatten()
    yramp1f = yramp1.flatten()
    reg_w = sp.interpolate.interp2d(yrampf,xrampf,Cwarp, yramp1f, xramp1f,'cubic');


def warpRegion(imgRegion, roi_pad, cx, cy, alpha):
    #Unpack region indices
    sx, ex, sy, ey = roi_pad

    x_nowarp, y_nowarp = np.mgrid[sx:ex+1, sy:ey+1]

    x_warped = alpha*(x_nowarp - cx)+cx;
    y_warped = alpha*(y_nowarp - cy)+cy;
    
    warped_image=
    # yrampc = yramp - cy;
    # # xramp1 = 1/wz*xrampc;
    # # yramp1 = 1/wz*yrampc;
    # xramp1 =alpha*xrampc;
    # yramp1 = 1/wz*yrampc;
    xrampf = xrampc.flatten()
    yrampf = yrampc.flatten()
    xramp1f = xramp1.flatten()
    yramp1f = yramp1.flatten()
    reg_w = sp.interpolate.interp2d(x_warped.flatten(),y_warped.flatten(),Cwarp, yramp1f, xramp1f,'cubic');



ldimage()
analyzeCA('full', im)