#!/usr/bin/env python

"""
2015 Jun 5 -- Runs the K2 SuPerPiG search
  Started by Brian Jackson decaelus@gmail.com
  Modified by Elisabeth Adams adams@psi.edu
"""
import numpy as np
#import matplotlib.pyplot as pl
import pylab as pl

#the fortran version of eebls
from eebls import eebls
#the MUCH slower python version -- don't use it!
#from eebls_python import eebls

import time as systime
import os.path
import astropy.io.fits as fits


def main(args):
    """
    Do the search
    """
    fits_files, aper, klobber, recondition, highersnr = args

    #Generate EPIC number list from fits
    epic_num_list = list()

    print("#epic_num, bper, bpow, depth, qtran, in1, in2")

    for current_file in fits_files:
        
        #Remove whitespace
        current_file = current_file.strip()
       
        #load in file
        hdulist = fits.open(current_file)
        #BESTAPER is the photometry from Vanderburg's best-fit aperture
        #  -- https://archive.stsci.edu/prepds/k2sff/
        HDU = hdulist[aper]

        #Get object EPIC number 
        epic_num = get_epic_num(HDU=HDU)

        epic_num_list.append(epic_num)

        # Improve work flow (2015 Aug 7 era) to check if conditioning has already been done
        if recondition == False:
            conditioned_data_file_name = make_output_file_name(current_file, replacement="conditioned_data.txt")
            detrended_data = np.loadtxt(conditioned_data_file_name, delimiter=' ')
            time = detrended_data[:,0]
            detrended_flux = detrended_data[:,1]
   
        else:
            #Load time series
            time, flux = get_time_series(HDU)

            #apply median boxcar filter to time series
            if highersnr: ### For transits we identified as very deep, we don't want to cutoff the bottom
                time, detrended_flux = detrend_flux(time, flux, outlier_threshold=100.)
            else:
                time, detrended_flux = detrend_flux(time, flux)
           # print min(detrended_flux),  np.median(detrended_flux), max(detrended_flux)

            #2015 Jun 26 -- Mask out thruster firing signal using periodic signal
            #time, detrended_flux = mask_thruster_events(time, detrended_flux)
            #2015 Sep 2 -- mask it directly from the MOVING flag in the fits file -- no need for additional function call

            #save detrended data to file
            conditioned_data_file_name = \
                write_conditioned_data_file(current_file, time, detrended_flux)


        # At this point the data has most definitely been conditioned!
        #apply eebls
        freqs, eebls_spec, bper, bpow, depth, qtran, in1, in2 = \
                apply_eebls(time, detrended_flux)

        #write out results to screen -- should be piped into results file on command line
        #Remember that Fortran starts indexing arrays at 1, while python starts at 0!!
        print("%s, %1.10f, %1.10f, %1.10f, %f, %i, %i" % (epic_num, bper, bpow, depth, qtran, in1-1, in2-1))

        #write eebls spectrum file
        eebls_spec_file_name = \
                write_eebls_spec_file(current_file, freqs, eebls_spec)

        #write spectrum plot to file
        eebls_spec_plot_file_name = \
                plot_eebls_spec(current_file, freqs, eebls_spec, bper, depth)

        # sometimes eebls fails, so don't try to fold!
        #make binned transit
        if bper == 0.0:
            binned_time, binned_flux, binned_err = \
                    make_binned_transit(time, detrended_flux)
            
        else:
            binned_time, binned_flux, binned_err = \
                    make_binned_transit(time % bper, detrended_flux)

        binned_time_double, binned_flux_double, binned_err_double = \
                    make_binned_transit(time % (2*bper), detrended_flux)
    

        #plot folded transit
        folded_transit_file_name = \
                plot_folded_transit(current_file, binned_time, binned_flux, binned_err,
                                    binned_time_double, binned_flux_double, binned_err_double, in1, in2, depth)

def get_now():
    now = systime.strftime("%Y%b%d")

    return now

# Old version that relies on periodicity (and generous cuts, removes ~17% of data)
def mask_thruster_events_from_periodicity(time, detrended_flux, 
             thruster_period=5.8845995/24.,
             t0=0.091379783494442951, t1=0.14249920382272482):
    #2015 Jun 26 -- Masks out data points that fall within thruster events

    folded_time = time % thruster_period
    ind = (folded_time < t0) | (folded_time > t1)

    return time[ind], detrended_flux[ind]

## New version that relies on MOVING tag (removes ~8% of data)
#def mask_thruster_events(HDU, time, flux):
#    """
#    Remove points have been flagged as 'MOVING'
#    """
#
#    moving_flag = HDU.data.field('MOVING') # 1 = yes it was moving (don't use!)
#    mask = np.array(1-np.array(moving_flag), dtype=bool) # sign flip: 1 = yes use point
#    time_masked = np.array(time)[mask]
#    flux_masked = np.array(flux)[mask]
#
#    return time_masked, flux_masked
#    

def plot_folded_transit(fits_file, binned_time, binned_flux, binned_err, 
                        binned_time_double, binned_flux_double, binned_err_double,
                        in1, in2, depth):

    now = get_now()

    fig = pl.figure(figsize=(5,5))
    ax = fig.add_subplot(211)
    ax.scatter(binned_time, binned_flux, marker='.')
    ax.set_ylim([min(binned_flux)-2*depth, max(binned_flux)+2*depth])
    ax.set_xlim([min(binned_time), max(binned_time)])
    ax.plot([min(binned_time), max(binned_time)],[0,0],ls="--",color="k")
    ax.plot([min(binned_time), max(binned_time)],[-depth,-depth],ls="--",color="k")
    pl.title("EPIC" + get_epic_num(filename=fits_file))
    pl.ylabel("$\Delta\ F$", fontsize=24)
    
    ax = fig.add_subplot(212)
    ax.scatter(binned_time_double, binned_flux_double, marker='.')
    ax.set_ylim([min(binned_flux)-depth, max(binned_flux)+depth])
    ax.set_xlim([min(binned_time_double), max(binned_time_double)])
    ax.plot([min(binned_time_double), max(binned_time_double)],[0,0],ls="--",color="k")
    ax.plot([min(binned_time_double), max(binned_time_double)],[-depth,-depth],ls="--",color="k")
    pl.xlabel("$t$ (days)", fontsize=24)
 
    binned_transit_plot_file_name = make_output_file_name(fits_file, 
        replacement="binned_transit_" + now + ".png")

    pl.savefig(binned_transit_plot_file_name, dpi=250, bbox_inches="tight")
    pl.close(fig)

    return binned_transit_plot_file_name

def make_binned_transit(time, detrended_flux, nb=200):
    import scipy.stats as ss
    from statsmodels.robust.scale import mad

    bins = np.linspace(min(time), max(time), num=nb)
    binned_time = np.array([0.5*(bins[i+1] + bins[i]) 
        for i in range(len(bins)-1)])
    binned_flux = ss.binned_statistic(time, detrended_flux, bins=bins)[0]
    binned_err = 1.4826*ss.binned_statistic(time, detrended_flux, bins=bins, 
            statistic=mad)[0]

    return binned_time, binned_flux, binned_err

def make_output_file_name(fits_file, replacement="eebls_spec.png"):
    return fits_file.replace("llc.fits", replacement)

def write_conditioned_data_file(fits_file, time, detrended_flux):
    now = get_now()

    conditioned_data_file_name = make_output_file_name(fits_file, 
            replacement="conditioned_data_" + now + ".txt")

    #write spectrum to file
    np.savetxt(conditioned_data_file_name, zip(time, detrended_flux),
            header="EPIC%s\n%s\n%s\n" % 
            (get_epic_num(filename=fits_file), 
                (fits_file.split('/'))[6], now), fmt="%f %1.10f") # was truncating too early!

    return conditioned_data_file_name

def plot_eebls_spec(fits_file, freqs, eebls_spec, bper, depth):

    now = get_now()

    fig = pl.figure(figsize=(5,5))
    ax1 = fig.add_subplot(111)
    ind = eebls_spec > 0.
    ax1.plot(24./freqs[ind], eebls_spec[ind], color='blue')

    ax1.axvline(bper*24., color='red', ls='--', lw=6)
    ax1.plot(24./freqs[ind], eebls_spec[ind], color='blue')

    pl.xlabel("$P$ (hours)", fontsize=24)
    pl.ylabel("BLS spec", fontsize=24)

    pl.xlim([min(24./freqs[ind]), max(24./freqs[ind])])
    pl.ylim([min(eebls_spec[ind]), max(eebls_spec[ind])])

    pl.title("EPIC" + get_epic_num(filename=fits_file))

    eebls_spec_plot_file_name = make_output_file_name(fits_file,
        replacement="eebls_spec_" + now + ".png")

    pl.savefig(eebls_spec_plot_file_name, dpi=250, bbox_inches="tight")

    pl.close(fig)

    return eebls_spec_plot_file_name

def write_eebls_spec_file(fits_file, freqs, eebls_spec):
    """
    Writes eebls spectrum to file
    """

    now = get_now()

    #make output file name
    output_file_name = make_output_file_name(fits_file,
            replacement="eebls_spec_" + now + ".txt")

    #write spectrum to file
    np.savetxt(output_file_name, zip(freqs, eebls_spec), 
            header="EPIC%s\n%s\n%s\n" % 
            (get_epic_num(filename=fits_file), (fits_file.split('/'))[6], now),
            fmt="%f %1.10f") # again need precision!

    return output_file_name

def apply_eebls(time, flux, num_periods=10000, 
        min_period=3./24, max_period=3., nb=100, qmi=0.01, qma=0.5):
    # min_period was 3 hours and max_period was 1 day
    # qmi was 0.01 and qma was 0.5

    #Set up the variables that eebls wants
    t = time
    x = flux
    n = len(t)
    u = np.zeros_like(t)
    v = np.zeros_like(t)
    
    #Convert periods to frequencies
    fmin = 1./max_period
    fmax = 1./min_period
    nf = num_periods
    df = (fmax - fmin)/(num_periods - 1)
    freqs = np.linspace(fmin, fmax, nf)

    #the fortran version of eebls
    p,bper,bpow,depth,qtran,in1,in2 = eebls(t,x,u,v,nf,fmin,df,nb,qmi,qma,n)

 #   # Well, this probably isn't required now htat rolling windows are used for sigma calculations
#    if bper == 0.0:
 #       qmi = 2.0*qmi
 #       print "zero period; incresing qmi to",qmi
 #       p,bper,bpow,depth,qtran,in1,in2 = eebls(t,x,u,v,nf,fmin,df,nb,qmi,qma,n)
  

    #The version of the call to use if using eebls_python,
    #  which is ridiculously slower than the fortran version, so don't.
#   p,bper,bpow,depth,qtran,in1,in2 = eebls(n,t,x,u,v,nf,fmin,df,nb,qmi,qma)

    return freqs, p, bper, bpow, depth, qtran, in1, in2

def detrend_flux(time, flux, window=1., outlier_threshold=10.):
    """
    Apply median boxcar filter to time series

    window -- window in days over which to apply boxcar filter

    """

    from scipy.signal import medfilt

    #First calculate the sampling frequency
    sampling = np.nanmedian([time[i+1] - time[i] for i in range(len(time)-1)])
    window_int = int(np.floor(1./sampling))

    #Make sure window is odd
    if(window_int % 2 == 0):
        window_int += 1

    #Since median filter from scipy assumes regular time sampling,
    #  the first thing I have to do is to interpolate the data to a 
    #  regular sampling grid. Hopefully, this step doesn't distort the data
    #  significantly.

    #generate light curve linearly interpolated to regular time sampling
    interp_time = np.arange(time[0], time[-1], sampling)
    interp_flux = np.interp(interp_time, time, flux)

    #calculate median filter signal
    smooth_interp_signal = medfilt(interp_flux, kernel_size=window_int)

    #interpolate back to the original irregularly sampled light curve
    reinterp_smooth_signal = np.interp(time, interp_time, smooth_interp_signal)

    #And subtract out the median filter
    smoothed_flux = flux - reinterp_smooth_signal

    if(outlier_threshold is not None):
        time, smoothed_flux = drop_outliers(time, smoothed_flux, thresh=outlier_threshold)
        #time, smoothed_flux = drop_outliers_rolling(time, smoothed_flux, thresh=outlier_threshold)
        #time, smoothed_flux = drop_outliers_smooth(time, smoothed_flux, reinterp_smooth_signal, thresh=outlier_threshold)
 
    return time, smoothed_flux

def drop_outliers(time, flux, thresh=10.):
    med = np.nanmedian(flux)
    sig = robust_stddev(flux)

    ind = abs(flux - med)/sig < thresh

    return time[ind], flux[ind]

def drop_outliers_smooth(time, flux,smooth_curve, thresh=10.):
    med = np.median(flux)
    sig = robust_stddev(flux)

    ind = abs(flux - med)/sig < thresh

    return time[ind], flux[ind]


def drop_outliers_rolling(time, flux, rolling_window=10, thresh=10.):
    ind=[]
    for ii in range(len(flux)):
        if ii <= 2*rolling_window:
            window = list(np.arange(0, 2*rolling_window,1))
        elif (len(flux) - ii) <= 2*rolling_window:
            window = list(np.arange(len(flux)-2*rolling_window,len(flux),1))
        else:
            window = list(np.arange(ii-rolling_window, ii+rolling_window,1))

        med = np.nanmedian(flux[window])
        sig = robust_stddev(flux[window])
        ind.append( abs(flux[ii] - med)/sig < thresh )
    ind = np.array(ind)

    return time[ind], flux[ind]


def robust_stddev(data):
    from statsmodels.robust.scale import mad as med_abs_dev

    #http://en.wikipedia.org/wiki/Median_absolute_deviation
    return 1.4826*med_abs_dev(data)

def get_time_series(HDU, which_flux='FCOR',mask_flag='MOVING',ditch_bad=True):
    """
    Get time series from data unit
    Remove points with mask_flag='MOVING' (or other keyword; None turns off masking)
    """

    #Add K2 time offset
    t0 = HDU.header['BJDREFI']

    time = HDU.data.field('T') + t0
    flux = HDU.data.field(which_flux)

    if mask_flag == None: # allow masking to be turned off if desired
        time_masked = time
        flux_masked = flux
    else:
        moving_flag = HDU.data.field(mask_flag) # 1 = yes it was moving (don't use!)
        mask = np.array(1-np.array(moving_flag), dtype=bool) # sign flip: 1 = yes use point
        time_masked = np.array(time)[mask]
        flux_masked = np.array(flux)[mask]
    # Only one case so far: -1 for most of 206166251 (c03)
    # Also two whole light curves: 202137146 and 202093417
    if ditch_bad:
        not_neg_ones = (flux_masked != -1.)
        time_masked = time_masked[not_neg_ones]
        flux_masked = flux_masked[not_neg_ones]
 
    return time_masked, flux_masked
    
 

def get_epic_num(HDU=None, filename=None):
    
    if(HDU is not None):
        return ((HDU.header['OBJECT']).split())[1]
    elif(filename is not None):
        ind1 = filename.index("lightcurve_")
        return filename[ind1 + 11: ind1 + 20]
    else:
        raise ValueError

def retrieve_files(dir):
    """
    Retrieves FITS files from root directory
    """
    import glob

    return glob.glob(dir + "*fits")

if __name__ == "__main__":
    
    import sys, argparse, glob

    #retrieve command-line params
    parser = argparse.ArgumentParser(description='Detrend K2 data and search for planetary transits using EEBLS.')
    parser.add_argument('--file', metavar='file', type=str, nargs=1, 
            help='single file on which to run detrend and eebls search')
    parser.add_argument('--fl', metavar='fl', type=str, nargs=1, 
            help='Name of a text file containing a list of all the fits files to process')
    parser.add_argument('--dir', metavar='dir', type=str, nargs=1, 
            default='/Users/bjackson/research/k2sff/c02/',
            help='root directory in which to find FITS data files')
    parser.add_argument('--aper', metavar='aper', type=str, nargs=1, 
            default='BESTAPER', help='which of Vanderburg\'s aperatures to use')

    parser.add_argument('--recondition', dest='recondition', action='store_true')
    parser.add_argument('--highersnr', dest='highersnr', action='store_true')

    parser.add_argument('--klobber', dest='klobber', action='store_true')
    parser.add_argument('--no-klobber', dest='klobber', action='store_false')
    parser.set_defaults(feature=False, recondition=True)

    args = parser.parse_args()

    single_file = args.file
    file_list = None
    if(args.fl is not None):
        file_list = args.fl[0]
    aper = args.aper
    klobber = args.klobber
    highersnr = args.highersnr
    recondition = args.recondition
    root_dir = args.dir

    if((file_list is not None) and (single_file is None)):
        f = open(file_list, 'r')
        files = f.readlines()
        f.close()
    elif((file_list is None) and (single_file is not None)):
        files = single_file
    else:
        files = glob.glob(root_dir + "*.fits")

    main([files, aper, klobber, recondition, highersnr])
