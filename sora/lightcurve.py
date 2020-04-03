import numpy as np
import matplotlib.pylab as pl
import astropy.units as u
from astropy.timeseries import BoxLeastSquares
from astropy.time import Time
import scipy.special as scsp
from datetime import datetime
from .extra import ChiSquare
import os

def calc_fresnel(distance,lambida):
    """ Returns the fresnel scale.
    ----------
    Parameters
    ----------
    distance  (int, float, array): distances, in km.
    lambida   (int, float, array): wavelength, in km.
    ----------
    """
    return np.sqrt(lambida*distance/2)


class LightCurve():
    '''
    Docstring
    Define a Light Curve
    '''
    def __init__(self, **kwargs):
        input_done = False
        self.dflux = None
        if 'file' in kwargs:
            if not os.path.isfile(kwargs['file']):
                raise ValueError('{} not found'.format(kwargs['file']))
            try:
                time, self.flux, self.dflux = np.loadtxt(kwargs['file'], usecols=[0,1,2], unpack=True)
            except:
                time, self.flux = np.loadtxt(kwargs['file'], usecols=[0,1], unpack=True)
            else:
                raise ValueError('Input file must ')
            input_done = True
        if 'time' in kwargs and 'flux' in kwargs:
            if input_done:
                raise ValueError('Only one type of input can be given. Please refer to the tutorial.')
            self.flux = kwargs['flux']
            time = kwargs['time']
            if len(self.flux) != len(time):
                raise ValueError('time and flux must have the same length')
            if 'dflux' in kwargs:
                self.dflux = kwargs['dflux']
                if len(self.flux) != len(self.dflux):
                    raise ValueError('dflux must have the same length as flux and time')
            input_done = True
        if input_done and 'exptime' not in kwargs:
            raise ValueError('exptime not defined')
        if 'exptime' in kwargs:
            self.exptime = kwargs['exptime']
        if 'tref' in kwargs:
            try:
                if type(kwargs['tref']) in [Time,str]:
                    self.tref = Time(kwargs['tref'])
                elif type(kwargs['tref']) in [int,float]:
                    self.tref = Time(kwargs['tref'], format='jd')
            except:
                raise ValueError('tref must be an Julian Date or ISO format Date')
        if 'time' in locals():
            if type(time) == Time:
                if not 'tref' in kwargs:
                    self.tref = Time(time[0].iso.split(' ')[0] + ' 00:00:00.000')
            elif all(time > 2400000):
                time = Time(time, format='jd')
                if not 'tref' in kwargs:
                    self.tref = Time(time[0].iso.split(' ')[0] + ' 00:00:00.000')
            elif not hasattr(self,'tref'):
                raise ValueError('tref must be given')
            else:
                time = self.tref + time*u.s
            self.time = (time - self.tref).sec
            order = np.argsort(self.time)
            self.flux = self.flux[order]
            self.time = self.time[order]
            if self.dflux != None:
                self.dflux = self.dflux[order]
            self.initial_time = np.min(time)
            self.end_time = np.max(time)
            self.cycle = np.median(self.time[1:] - self.time[:-1])
        if 'immersion' in kwargs and 'emersion' in kwargs:
            self._immersion = kwargs['immersion']
            self._emersion = kwargs['emersion']
            self.immersion_err = 0.0
            if 'immersion_err' in kwargs:
                self.immersion_err = kwargs['immersion_err']
            self.emersion_err = 0.0
            if 'emersion_err' in kwargs:
                self.emersion_err = kwargs['emersion_err']
            input_done = True
        if 'initial_time' in kwargs and 'end_time' in kwargs:
            self.initial_time = kwargs['initial_time']
            self.end_time = kwargs['end_time']
            input_done = True
        if not input_done:
            raise ValueError('No allowed input conditions satisfied. Please refer to the tutorial.')
        self.lambda_0 = 0.70 #microns #### *u.micrometer.to('km')
        if 'lambda' in kwargs:
            self.lambda_0 = kwargs['lambda']
        self.delta_lambda = 0.30 #microns #### *u.micrometer.to('km')
        if 'delta_lambda' in kwargs:
            self.delta_lambda = kwargs['delta_lambda']
        if 'dist' in kwargs:
            self.set_dist(kwargs['dist'])
        if 'diam' in kwargs:
            self.set_diam(kwargs['diam'])
        if 'vel' in kwargs:
            self.set_vel(kwargs['vel'])
        if hasattr(self,'time'):
            self.model = np.ones(len(self.time))
        self.dt = 0.0

    @property
    def immersion(self):
        return self._immersion + self.dt*u.s

    @property
    def emersion(self):
        return self._emersion + self.dt*u.s
    
    def set_vel(self,vel):
        '''
        Set the occultation velocity
        Inputs:
        vel = float, in km/s
        '''
        if type(vel) == u.quantity.Quantity:
            vel = vel.to(u.km/u.s).value
        elif type(vel) == [float,int]:
            pass
        else:
            raise TypeError('vel must be a float or an Astropy Unit object')
        self.vel = np.absolute(vel)

    def set_dist(self,dist):
        '''
        Set the object distance
        Inputs:
        dist = float, in km
        '''
        if type(dist) == u.quantity.Quantity:
            dist = dist.to(u.AU).value
        elif type(dist) in [float,int]:
            pass
        else:
            raise TypeError('dist must be a float or an Astropy Unit object')
        self.dist = dist

    def set_diam(self,diam):
        '''
        Set the star diameter
        Inputs:
        diam = float, in km
        '''
        if type(diam) == u.quantity.Quantity:
            diam = diam.to(u.km).value
        elif type(diam) in [float,int]:
            pass
        else:
            raise TypeError('diam must be a float or an Astropy Unit object')
        self.d_star = diam
    
    def occ_model(self,t_ingress, t_egress, opa_ampli, mask, npt_star=12, time_resolution_factor=10):
        """ Returns the modelled light curve considering fresnel difraction, star diameter and intrumental response.
        ----------
        Parameters
        ----------
        t_ingress (int, float): Ingrees time, in seconds.                       (input)
        t_egress  (int, float): Egress time, in seconds.                        (input)
        opa_ampli (int, float): Opacity, opaque = 1.0, transparent = 0.0.       (input)
        mask (array with Booleans): Mask with True values to be computed        (input)
        npt_star  (int): Number of subdivisions for computing the star size's effects, default equal to 12. (auto)
        time_resolution_factor (int,float): Steps for fresnel scale used for modelling the light curve,     (auto)
            default equals to 10 steps for fresnel scale.
        ----------
        Returns
        ----------
        self.flux_inst       (array): Modelled Instrumental light flux.
        self.time_model      (array): Modelled timing.
        self.flux_star       (array): Modelled light flux considering fresnel difraction and star's diameter.
        self.flux_fresnel    (array): Modelled light flux considering fresnel difraction.
        self.model_geometric (array): Modelled light flux considering a box model.
        """
        # Computing the fresnel scale
        lamb  = self.lambda_0*u.micrometer.to('km')
        dlamb = self.delta_lambda*u.micrometer.to('km')
        dist  = self.dist*u.au.to('km')
        vel = np.absolute(self.vel)
        time_obs = self.time[mask]
        fresnel_scale_1 = calc_fresnel(dist,lamb-dlamb/2.0)
        fresnel_scale_2 = calc_fresnel(dist,lamb+dlamb/2.0)
        fresnel_scale   = (fresnel_scale_1 + fresnel_scale_2)/2.0
        time_resolution = (np.min([fresnel_scale/vel,self.exptime]))/time_resolution_factor
        #
        #Creating a high resolution curve to compute fresnel difraction, stellar diameter and instrumental integration 
        time_model = np.arange(time_obs.min()-5*self.exptime,time_obs.max()+5*self.exptime,time_resolution)
        #
        #Changing X: time (s) to distances in the sky plane (km), considering the tangential velocity (vel in km/s)
        x   = time_model*vel    
        x01 = t_ingress*vel
        x02 = t_egress*vel
        #
        #Computing fresnel diffraction for the case where the star size is negligenciable
        flux_fresnel_1 = self.__bar_fresnel(x,x01,x02,fresnel_scale_1,opa_ampli)
        flux_fresnel_2 = self.__bar_fresnel(x,x01,x02,fresnel_scale_2,opa_ampli)
        flux_fresnel   = (flux_fresnel_1 + flux_fresnel_2)/2.
        flux_star      = flux_fresnel.copy()
        if (self.d_star > 0):
            #Computing fresnel diffraction for the case where the star size is not negligenciable
            resolucao   = self.d_star/npt_star
            flux_star_1 = np.zeros(len(time_model))
            flux_star_2 = np.zeros(len(time_model))
            #Computing stellar diameter only near the ingrees or egrees times
            star_diam = (np.absolute(x - x01) < 3*self.d_star) + (np.absolute(x - x02) < 3*self.d_star)
            p = np.arange(-npt_star,npt_star)*resolucao
            coeff = np.sqrt(np.absolute(self.d_star**2 - p**2))
            for ii in np.where(star_diam == True)[0]:
                xx = x[ii] + p
                flux1 = self.__bar_fresnel(xx,x01,x02,fresnel_scale_1,opa_ampli)
                flux2 = self.__bar_fresnel(xx,x01,x02,fresnel_scale_2,opa_ampli)
                flux_star_1[ii] = np.sum(coeff*flux1)/coeff.sum()
                flux_star_2[ii] = np.sum(coeff*flux2)/coeff.sum()
                flux_star[ii]   = (flux_star_1[ii] + flux_star_2[ii])/2.
        flux_inst = np.zeros(len(time_obs)) 
        for i in range(len(time_obs)):
            event_model  = (time_model > time_obs[i]-self.exptime/2.) & (time_model < time_obs[i]+self.exptime/2.)
            flux_inst[i] = (flux_star[event_model]).mean()
        self.model[mask] = flux_inst
        self.time_model = time_model
        self.model_star = flux_star
        self.model_fresnel = flux_fresnel
        ev_model = (time_model > t_ingress) & (time_model < t_egress)
        flux_box = np.ones(len(time_model))
        flux_box[ev_model] = (1-opa_ampli)**2
        self.model_geometric = flux_box
        return


  
    #def occ_lcfit(self,mask, t_ingress, t_egress, opa_ampli=1, dt_ingress=0, dt_egress=0, dopacity=0, loop=10000):
    def occ_lcfit(self, **kwargs):
        """ Brute force chi square fit for occultations lightcurve.
        ----------
        Parameters
        ----------
        tmin (int,float): Minimum time to consider in the fit prcedure                  (input)
        tmax (int,float): Maximum time to consider in the fit prcedure                  (input)
        t_ingress  (int, float): Ingrees time, in seconds.                              (input)
        t_egress   (int, float): Egress time, in seconds.                               (input)
        opacity  (int, float): Opacity, opaque = 1.0, transparent = 0.0, default = 1. (input)
        delta_t (int, float): Interval to fit ingress or egress time                       (auto)
        dopacity   (int, float): Interval to fit opacity, default equal to 0, no fit.      (auto)
        loop       (int): Number of tests to be done, default equal to 10000.              (auto)
        ----------
        Returns
        ----------
        chi2 (ChiSquare): ChiSquare object
        """
        #
        if not hasattr(self, 'flux'):
            raise ValueError('Fit curve is only possible when a LightCurve is instatiated with time and flux.')
        delta_t = 2*self.cycle
        loop = 10000
        t_i = np.zeros(loop)
        t_e = np.zeros(loop)
        tmax = self.time.max()
        tmin = self.time.min()
        t_ingress = tmin - self.exptime
        do_ingress = False
        t_egress = tmax + self.exptime
        do_egress = False
        opacity = 1.0
        delta_opacity = 0.0
        do_opacity = False
        if 'loop' in kwargs:
            loop = kwargs['loop']
        if 't_ingress' not in kwargs and 't_egress' not in kwargs:
            preliminar_occ = self.occ_detect()
            t_ingress = preliminar_occ['imersion_time']
            do_ingress = True
            t_egress = preliminar_occ['emersion_time']
            do_egress = True
            delta_t = 5*preliminar_occ['time_err']
            tmax = t_egress+2*preliminar_occ['occultation_duration']
            tmin = t_ingress-2*preliminar_occ['occultation_duration']
            if 2*preliminar_occ['occultation_duration'] < 10*self.cycle:
                tmax = t_egress + 10*self.cycle
                tmin = t_ingress - 10*self.cycle
        if 'tmax' in kwargs:
            tmax = kwargs['tmax']
        if 'tmin' in kwargs:
            tmin = kwargs['tmin']
        if 'delta_t' in kwargs:
            delta_t = kwargs['delta_t']
        if 't_ingress' in kwargs:
            t_ingress = kwargs['t_ingress']
            do_ingress = True
        t_i = t_ingress + delta_t*(2*np.random.random(loop) - 1)
        if 't_egress' in kwargs:
            t_egress = kwargs['t_egress']
            do_egress = True
        t_e = t_egress  + delta_t*(2*np.random.random(loop) - 1)
        mask = (self.time >= tmin) & (self.time <= tmax)
        mask_sigma = (((self.time >= tmin) & (self.time < t_ingress - self.exptime)) +
                      ((self.time > t_egress + self.exptime) & (self.time <= tmax)))
        sigma = self.flux[mask_sigma].std()
        if 'sigma' in kwargs:
            sigma = kwargs['sigma']
        if 'opacity' in kwargs:
            opacity = kwargs['opacity']
        if 'dopacity' in kwargs:
            delta_opacity = kwargs['dopacity']
            do_opacity = True
        opas = opacity + delta_opacity*(2*np.random.random(loop) - 1)
        opas[opas>1.], opas[opas<0.] = 1.0, 0.0
        #
        tflag = np.zeros(loop)
        tflag[t_i > t_e] = t_i[t_i > t_e]
        t_i[t_i > t_e]    = t_e[t_i > t_e]
        t_e[t_i > t_e]    = tflag[t_i > t_e]
        chi2 = 999999*np.ones(loop)
        #tcontrol_f0 = datetime.now()
        for i in range(loop):
            model_test = self.__occ_model(t_i[i],t_e[i],opas[i],mask)
            chi2[i] = np.sum((self.flux[mask] -  model_test)**2)/(sigma**2)
        #tcontrol_f3 = datetime.now()
        #print('Elapsed time: {:.3f} seconds.'.format((tcontrol_f3 - tcontrol_f0).total_seconds()))
        kkargs = {}
        if do_ingress:
            kkargs['immersion'] = t_i
        if do_egress:
            kkargs['emersion'] = t_e
        if do_opacity:
            kkargs['opacity'] = opas
        chisquare = ChiSquare(chi2, len(self.flux[mask]), **kkargs)
        onesigma = chisquare.get_nsigma(1)
        if 'immersion' in onesigma:
            self._immersion = self.tref + onesigma['immersion'][0]*u.s
            self.immersion_err = onesigma['immersion'][1]
            t_ingress = onesigma['immersion'][0]
        if 'emersion' in onesigma:
            self._emersion = self.tref + onesigma['emersion'][0]*u.s
            self.emersion_err = onesigma['emersion'][1]
            t_egress = onesigma['emersion'][0]
        if 'opacity' in onesigma:
            opacity = onesigma['opacity'][0]
        # Run occ_model() to save best parameters in the Object.
        self.occ_model(t_ingress,t_egress,opacity,mask)
        return chisquare

    def plot_lc(self):
        '''
        Plot the light curve if you want to save the file, the fig_name should be different than None
        '''
        if (np.any(self.time) != None):
            pl.close()
            pl.figure(figsize=[8.4, 3.4])
            pl.plot(self.time,self.flux,'k.-',label='Obs.',zorder=1)
            if (np.any(self.model) != None):
                pl.plot(self.time,self.model,'r-',label='Model',zorder=2)
                pl.scatter(self.time,self.model, s=50, facecolors='none', edgecolors='r',zorder=3)
            pl.tight_layout()
            pl.xlabel('Time [seconds]',fontsize=20)
            pl.ylabel('Relative Flux',fontsize=20)
            pl.legend(fontsize=20,ncol=2)
            pl.xticks(fontsize=20)
            pl.yticks(fontsize=20)
        else:
            raise ValueError('There is no lightcurve to plot')
        return
    
    def occ_detect(self, maximum_duration=None, dur_step=None, snr_limit=None, \
                  n_detections=None):
        """
        Detect automatically the occultation event in the light curve
        
        Parameters
        ----------
        maximum_duration: float, optional
            Maximum duration of the occultation event (default is 1/4th of the
            light curve's time span).
        dur_step: float, optionl
            Step size to sweep occultation duration event (default value is 1/2
            of sampling).
        snr_limit: float, optional
            Minimum occultation SNR.
        n_detections: int, optional
            N best detections regardless from SNR. n_detections is superseded
            by snr_limit.
            
        Returns
        -------
        OrderedDict
            An ordered dictionary of :attr:`name`::attr:`value` pairs for each
            Parameter.
            
        Examples
        --------
        >>> params = lc.occ_detect()
        >>> params
        {'occultation_duration': 0.0004645648878067732,
         'central_time': 2457852.5916293273,
         'imersion_time': 2457852.5913970447,
         'emersion_time': 2457852.59186161,
         'time_err': 5.799811333417892e-07,
         'depth': 0.8663887801707082,
         'depth_err': 0.10972550419008305,
         'baseline': 0.9110181732552853,
         'baseline_err': 0.1904360360568157,
         'snr': 91.21719495827487,
         'occ_mask': array([False, False, False, ..., False, False, False])}
        """

        # duration of the light curve
        time_span = self.time[-1]-self.time[0]
        if maximum_duration and (maximum_duration > time_span):
            raise ValueError('Occultation duration (maximum_duration={0}) ' \
                             'exceeds the time series lenght ({1:0.5f}).' \
                             .format(maximum_duration, time_span))
        if not maximum_duration:
            maximum_duration = time_span*0.25
        if not dur_step:
            dur_step = np.median(self.time[1:-1]-self.time[0:-2])/2      
        duration_grid = np.arange(dur_step, maximum_duration, dur_step)
        # initial occultation mask (all data points)
        mask = np.ones(len(self.time), dtype=bool)

        if snr_limit:
            # minimum SNR accepted in a detection for multiple search
            snr_value = snr_limit+1
            occ0 = self.__run_bls(time_span, duration_grid)
            mask *= ~occ0['occ_mask']
            while (snr_value > snr_limit):
                occ1 = self.__run_bls(time_span, duration_grid, mask=mask)
                if occ1['snr'] > snr_limit:
                    snr_value = occ1['snr']
                    mask *= ~occ1['occ_mask']                    
                    occ0 = self.__summarize_bls(occ0,occ1)
                else:
                    snr_value = snr_limit
            return occ0
        elif n_detections:
            # search the n best fits
            occ0 = self.__run_bls(time_span, duration_grid)
            mask *= ~occ0['occ_mask']    
            for i in range(n_detections-1):
                occ1 = self.__run_bls(time_span, duration_grid, mask=mask)
                snr_value = occ1['snr']
                mask *= ~occ1['occ_mask']                    
                occ0 = self.__summarize_bls(occ0,occ1)
            return occ0
        else:
            # search only the first best fit
            return self.__run_bls(time_span, duration_grid)
        
        
    def __run_bls(self, per_grid, dur_grid, mask=None):
        """
        Private function to find the best box fit suitable to the data
        """

        # object with no occultation mask
        mskmodel = BoxLeastSquares(self.time, self.flux, dy=self.dflux)
        # if there is no dflux array, reset it to None in case of 
        # using a mask
        if self.dflux is None:
            dfluxmask = None
        else:
            dfluxmask = self.dflux[mask]
            
        # object with occultation mask
        if np.sum(mask):
            model = BoxLeastSquares(self.time[mask], self.flux[mask], \
                                    dy=dfluxmask)
        else:
            model = mskmodel   
        
        r = model.power(per_grid, dur_grid, objective='snr', method='fast')
        # statistics of the BLS fit
        stats = model.compute_stats(r.period, r.duration, r.transit_time)
        # occultation mask of the event with respect to all data
        occ_mask = mskmodel.transit_mask(self.time, r.period, r.duration, \
                                         r.transit_time)
        # parameters computation for clarity purposes
        occultation_duration = r.duration[0]
        central_time = stats['transit_times'][0]
        imersion_time = stats['transit_times'][0] - r.duration[0]/2
        emersion_time = stats['transit_times'][0] + r.duration[0]/2
        time_err = dur_grid[1] - dur_grid[0]
        depth = np.mean(self.flux[~occ_mask])-np.mean(self.flux[occ_mask])
        depth_err = np.std(self.flux[occ_mask])
        baseline = np.mean(self.flux[~occ_mask])
        baseline_err = np.std(self.flux[~occ_mask])
        snr = (depth/baseline_err)*np.sqrt(np.sum(occ_mask))

        return {'occultation_duration' : occultation_duration, 
                     'central_time' : central_time, 
                     'imersion_time' : imersion_time, 
                     'emersion_time' : emersion_time,
                     'time_err' : time_err,
                     'depth' : depth, 'depth_err' : depth_err,
                     'baseline' : baseline, 'baseline_err' : baseline_err,
                     'snr' : snr,'occ_mask' : occ_mask}
    
    def __summarize_bls(self, dict1, dict2):
        ''' Private function to merge dictionaries returned by BLS and 
            keep values of common keys in list.
        '''
        dict3 = {}
        for key, value in dict1.items():
            if key == 'occ_mask':
                sz = int(np.size(dict1[key])/np.size(dict2[key]))
                if sz > 1:
                    l = [None]*(sz+1)
                    for i in range(sz):
                        l[i] = dict1[key][i]
                    l[i+1] = dict2[key]
                    dict3[key] = l
                else:
                    dict3[key] = [dict1[key],dict2[key]]
            else:
                dict3[key] = np.append(dict1[key],dict2[key])
        return dict3

    def __bar_fresnel(self,X,X01,X02,fresnel_scale,opa_ampli):
        """ Returns the modelled light curve considering fresnel difraction.
        ----------
        Parameters
        ----------
        X   (array): Array with time values converted in km using the event velocity.
        X01 (int, float): Ingrees time converted in km using the event velocity.
        X02 (int, float): Egress time converted in km using the event velocity.
        fresnel_scale (int, float): Fresnel scale.
        opa_ampli     (int, float): Opacity, opaque = 1.0, transparent = 0.0
        ----------
        Returns
        ----------
        flux_fresnel (array):
        """
        # Converting from km to units of fresnel scale
        x   = X/fresnel_scale
        x01 = X01/fresnel_scale
        x02 = X02/fresnel_scale
        # Fresnel difraction parameters 
        x1 = x - x01
        x2 = x - x02
        s1,c1 = scsp.fresnel(x1)
        s2,c2 = scsp.fresnel(x2)
        cc = c1 - c2
        ss = s1 - s2
        r_ampli = - (cc+ss)*(opa_ampli/2.)
        i_ampli =   (cc-ss)*(opa_ampli/2.)
        # Determining the flux considering fresnel difraction
        flux_fresnel = (1.0 + r_ampli)**2 + (i_ampli)**2
        return flux_fresnel


    def __occ_model(self,t_ingress, t_egress, opa_ampli, mask, npt_star=12, time_resolution_factor=10):
        """ Private function returns the modelled light curve considering fresnel difraction, star diameter and intrumental response, intended for fitting inside the self.occ_lcfit().
        ----------
        Parameters
        ----------
        t_ingress (int, float): Ingrees time, in seconds.                       (input)
        t_egress  (int, float): Egress time, in seconds.                        (input)
        opa_ampli (int, float): Opacity, opaque = 1.0, transparent = 0.0.       (input)
        mask (array with Booleans): Mask with True values to be computed        (input)
        npt_star  (int): Number of subdivisions for computing the star size's effects, default equal to 12. (auto)
        time_resolution_factor (int,float): Steps for fresnel scale used for modelling the light curve,     (auto)
            default equals to 10 steps for fresnel scale.
        ----------
        Returns
        ----------
        flux_inst       (array): Modelled Instrumental light flux.
        """
        # Computing the fresnel scale
        lamb  = self.lambda_0*u.micrometer.to('km')
        dlamb = self.delta_lambda*u.micrometer.to('km')
        dist  = self.dist*u.au.to('km')
        vel = np.absolute(self.vel)
        time_obs = self.time[mask]
        fresnel_scale_1 = calc_fresnel(dist,lamb-dlamb/2.0)
        fresnel_scale_2 = calc_fresnel(dist,lamb+dlamb/2.0)
        fresnel_scale   = (fresnel_scale_1 + fresnel_scale_2)/2.0
        time_resolution = (np.min([fresnel_scale/vel,self.exptime]))/time_resolution_factor
        #
        #Creating a high resolution curve to compute fresnel difraction, stellar diameter and instrumental integration 
        time_model = np.arange(time_obs.min()-5*self.exptime,time_obs.max()+5*self.exptime,time_resolution)
        #
        #Changing X: time (s) to distances in the sky plane (km), considering the tangential velocity (vel in km/s)
        x   = time_model*vel    
        x01 = t_ingress*vel
        x02 = t_egress*vel
        #
        #Computing fresnel diffraction for the case where the star size is negligenciable
        flux_fresnel_1 = self.__bar_fresnel(x,x01,x02,fresnel_scale_1,opa_ampli)
        flux_fresnel_2 = self.__bar_fresnel(x,x01,x02,fresnel_scale_2,opa_ampli)
        flux_fresnel   = (flux_fresnel_1 + flux_fresnel_2)/2.
        flux_star      = flux_fresnel.copy()
        if (self.d_star > 0):
            #Computing fresnel diffraction for the case where the star size is not negligenciable
            resolucao   = self.d_star/npt_star
            flux_star_1 = np.zeros(len(time_model))
            flux_star_2 = np.zeros(len(time_model))
            #Computing stellar diameter only near the ingrees or egrees times
            star_diam = (np.absolute(x - x01) < 3*self.d_star) + (np.absolute(x - x02) < 3*self.d_star)
            p = np.arange(-npt_star,npt_star)*resolucao
            coeff = np.sqrt(np.absolute(self.d_star**2 - p**2))
            for ii in np.where(star_diam == True)[0]:
                xx = x[ii] + p
                flux1 = self.__bar_fresnel(xx,x01,x02,fresnel_scale_1,opa_ampli)
                flux2 = self.__bar_fresnel(xx,x01,x02,fresnel_scale_2,opa_ampli)
                flux_star_1[ii] = np.sum(coeff*flux1)/coeff.sum()
                flux_star_2[ii] = np.sum(coeff*flux2)/coeff.sum()
                flux_star[ii]   = (flux_star_1[ii] + flux_star_2[ii])/2.
        flux_inst = np.zeros(len(time_obs)) 
        for i in range(len(time_obs)):
            event_model  = (time_model > time_obs[i]-self.exptime/2.) & (time_model < time_obs[i]+self.exptime/2.)
            flux_inst[i] = (flux_star[event_model]).mean()
        return flux_inst
