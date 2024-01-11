import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import scipy.stats
import warnings
import csv
from acoustipy.TMM import AcousticTMM
from acoustipy.Database import AcoustiBase


warnings.filterwarnings("ignore", category=RuntimeWarning)

class AcousticID():
    """
    Create a AcousticID object
    
    Description:
    ------------
    AcousticID is an optimization routine built on top of SciPy (https://scipy.org/) and AcousticTMM that can be used 
    to identify the difficult-to-measure parameters of the Johnson-Chompoux-Allard equivalent fluid model using an 
    inverse characterization procedure based on the following paper: 
    
    Atalla, Youssef & Panneton, R.. (2005). Inverse acoustical characterization of open cell porous 
    media using impedance tube measurements. Canadian Acoustics - Acoustique Canadienne. 33.

    Or by using an indirect characterization procedure based on the following paper:

    Panneton, R. & Salissou, Yacoubou. (2009). Indirect acoustical characterization of sound absorbing materials.. 
    The Journal of the Acoustical Society of America. 126. 2297. 10.1121/1.3249416.

    Or by a hybrid characterization procedure that combines the inverse and indirect methods.
    
    Attributes
    ----------
    mount_type (str):
        'No Gap', 'Gap', or 'Dual'    
        Specify whether impedance tube measurements are of a sample with rigid backing, an air gap, or both.
    
    no_gap_file (str):
        Name of the csv filepath that contains the frequency dependent absorption, reflection, or surface impedance
        coefficients of a single porous layer obtained from an impedence tube measurement with rigid backing.
        The csv file should contain 2 columns of equal length -- the frequencies in the 1st column and coefficients in the 2nd.
        
    gap_file (str):
        Name of the csv filepath that contains the frequency dependent absorption, reflection, or surface impedance
        coefficients of a single porous layer obtained from an impedence tube measurement with an air gap backing.
        The csv file should contain 2 columns of equal length -- the frequencies in the 1st column and coefficients in the 2nd.

    input_type (str):
        'absorption', 'reflection', or 'surface' -- specifies the type of measurement made with the impedance tube.

    air_temperature (float):
        Temperature of air [°C].  If specified, all other air properties will be determined by this parameter.
    
    sound_speed (float):
        Speed of sound in air [m/s]
    
    air_density (float):
        Density of air [kg/m3]
    
    Cp (float):
        Specific heat @ constant pressure [kJ/kg K]
    
    Cv (float):
        Specifc heat @ constant volume [kJ/kg K]
    
    viscosity (float):
        Dynamic viscosity of air [kg/m*s]
                                      
    Pr (float):
        Prandtl number of air []
    
    P0 (float):
        Atmospheric pressure [Pa]

    """
    def __init__(self,
                 mount_type:str = "No Gap",
                 no_gap_file:str = None,
                 gap_file:str = None,
                 input_type:str = 'absorption',
                 air_temperature:float = None,
                 sound_speed:float = 343.152,
                 air_density:float = 1.2058,
                 Cp:float = 1.004425,
                 Cv:float = 0.717425,
                 viscosity:float = 1.825e-05,
                 Pr:float = .7157,
                 P0:float = 101325
                 ):

        self.temp = air_temperature
        self.speed = sound_speed
        self.density = air_density
        self.Cp = Cp
        self.Cv = Cv
        self.viscosity = viscosity
        self.Pr = Pr
        self.P0 = P0
        self.opt_type = mount_type
        self.input_type = input_type
        self.thickness = None
        self.flow_resistivity = None
        self.porosity = None
        self.uncertainty = None
        self.air_gap = None
        
        

        if self.opt_type == 'No Gap':
            if self.input_type == 'reflection':
                self.no_gap_data = self.load_to_array(no_gap_file)
                          
            elif self.input_type == 'surface':
                self.no_gap_data = self.load_to_array(no_gap_file)
                
            elif self.input_type == 'absorption':
                    self.no_gap_data = self.load_to_array(no_gap_file,type='float')

        elif self.opt_type == 'Gap':
            if self.input_type == 'reflection':    
                self.gap_data = self.load_to_array(gap_file)
            
            elif self.input_type == 'surface':
                self.gap_data = self.load_to_array(gap_file)
                
            elif self.input_type == 'absorption':
                    self.gap_data = self.load_to_array(gap_file,type='float')

        elif self.opt_type == 'Dual':
            if self.input_type == 'reflection':
                self.no_gap_data = self.load_to_array(no_gap_file)
                self.gap_data = self.load_to_array(gap_file)
            
            elif self.input_type == 'surface':
                self.no_gap_data = self.load_to_array(no_gap_file)
                self.gap_data = self.load_to_array(gap_file)
                

            elif self.input_type == 'absorption':
                self.no_gap_data = self.load_to_array(no_gap_file,type='float')
                self.gap_data = self.load_to_array(gap_file,type='float')
        
    @property
    def frequency(self):
        #frequency range of interest
        
        if self.opt_type == 'No Gap':
            no_gap_freq = self.no_gap_data[:,0]
            return(no_gap_freq)
        
        elif self.opt_type == 'Gap':
            gap_freq = self.gap_data[:,0]
            return(gap_freq)
        
        elif self.opt_type == 'Dual':
            no_gap_freq = self.no_gap_data[:,0]
            gap_freq = self.gap_data[:,0]
        
            if np.array_equal(no_gap_freq,gap_freq) != True:
                raise ValueError("Frequencies must match between no gap and gap absorption curves")
            else:
                return(no_gap_freq)


    @property
    def meas_abs(self):
        #measured absorption coefficients
        if self.input_type == 'absorption':
            if self.opt_type == 'No Gap':
                no_gap_abs = self.no_gap_data[:,1]
                return([no_gap_abs,None])
            
            elif self.opt_type == 'Gap':
                gap_abs = self.gap_data[:,1]
                return([None,gap_abs])
            
            elif self.opt_type == 'Dual':
                no_gap_abs = self.no_gap_data[:,1]
                gap_abs = self.gap_data[:,1]
                return([no_gap_abs,gap_abs])

        elif self.input_type == 'reflection':
            if self.opt_type == 'No Gap':
                no_gap_abs = 1-abs(self.no_gap_data[:,1])**2
                return([no_gap_abs,None])
            
            elif self.opt_type == 'Gap':
                gap_abs = 1-abs(self.gap_data[:,1])**2
                return([None,gap_abs])
            
            elif self.opt_type == 'Dual':
                no_gap_abs = 1-abs(self.no_gap_data[:,1])**2
                gap_abs = 1-abs(self.gap_data[:,1])**2
                return([no_gap_abs,gap_abs])
        
        elif self.input_type == 'surface':
            if self.opt_type == 'No Gap':
                ngzs = self.no_gap_data
                ngr = (ngzs-self.Z0)/(ngzs+self.Z0)
                no_gap_abs = 1-abs(ngr[:,1])**2
                return([no_gap_abs,None])
            
            elif self.opt_type == 'Gap':
                gzs = self.gap_data
                gr = (gzs-self.Z0)/(gzs+self.Z0)
                gap_abs = 1-abs(gr[:,1])**2
                return([None,gap_abs])
            
            elif self.opt_type == 'Dual':
                ngzs = self.no_gap_data
                ngr = (ngzs-self.Z0)/(ngzs+self.Z0)
                no_gap_abs = 1-abs(ngr[:,1])**2

                gzs = self.gap_data
                gr = (gzs-self.Z0)/(gzs+self.Z0)
                gap_abs = 1-abs(gr[:,1])**2
                return([no_gap_abs,gap_abs])

    @property    
    def density_temp(self):
        #Temperature dependent density of air
        if self.temp is None:
            airdensity = self.density
        else:
            airdensity = ((1.07743*1e-5)*(self.temp**2))+((-0.004581339)*self.temp)+(1.294685259)
        return (airdensity)
    
    @property
    def soundspeed_temp(self):
        #Temperature dependent speed of sound in air
        if self.temp is None:
            soundspeed = self.speed
        else:
            soundspeed = (.6016*self.temp)+331.12
        return(soundspeed)
    
    @property
    def gamma_temp(self):
        #Temperature dependent specific heat ratio
        if self.temp is None:
            gamma = self.Cp/self.Cv
        else:
            Cp = (4.00166852057851e-07*self.temp**2)+(1.69769187986639e-05*self.temp)+(1.00559293937709)
            Cv = (3.65205412117683e-07*self.temp**2)+(2.88811127246258e-05*self.temp)+(7.17032243570935e-01)
            gamma = Cp/Cv
        return(gamma)
    
    @property
    def Cp_temp(self):
        #Temperature dependent specific heat ratio
        if self.temp is None:
            Cp = self.Cp
        else:
            Cp = (4.00166852057851e-07*self.temp**2)+(1.69769187986639e-05*self.temp)+(1.00559293937709)

        return(Cp)
    
    @property
    def Cv_temp(self):
        #Temperature dependent specific heat ratio
        if self.temp is None:
            Cv = self.Cv
        else:
            Cv = (3.65205412117683e-07*self.temp**2)+(2.88811127246258e-05*self.temp)+(7.17032243570935e-01)

        return(Cv)
    
    @property
    def viscosity_temp(self):
        #Temperature dependent viscosity of air
        if self.temp is None:
            visc = self.viscosity
            return(visc)
        else:
            visc = (-3.52159145081837e-11*self.temp**2)+(4.93272149610679e-8*self.temp)+(1.72293415521214e-5)
            return(visc)
    
    @property
    def Pr_temp(self):
        #Temperature dependent prandtl number
        if self.temp is None:
            Prandtl = self.Pr
        else:
            Prandtl = (7.06471476243448e-07*self.temp**2)+(-2.20826446051168e-04*self.temp)+(0.71980868)
        return(Prandtl)
    
    @property
    def Z0(self):
        #Characteristic Impedance of Air
        z0 = self.density_temp*self.soundspeed_temp
        return(z0)
        
    def _predictionJCA(self,
                       parameters: dict) -> np.ndarray:
        """
        Calculates the predicted frequency dependent absorption curve using the parameters identified in the
        optimization procedure defined in the find_values methods

        Parameters
        ----------
        parameters (dict):
            dictionary containing the identified thickness, flow resistivity, porosity, tortuosity, 
            viscous characteristic length, thermal characteristic length, and air gap of the sample.

        Returns
        -------
        predicted (ndarray):
            2D array of frequencies and predicted absorption coefficients 

        """
        
        dummy_struct = AcousticTMM(air_temperature=self.temp,
                                   sound_speed=self.soundspeed_temp,
                                   air_density=self.density_temp,
                                   Cp=self.Cp,
                                   Cv=self.Cv,
                                   viscosity=self.viscosity_temp,
                                   Pr=self.Pr_temp,
                                   P0=self.P0)

        t = parameters['thickness']
        fr = parameters['flow resistivity']
        phi = parameters['porosity']
        tort = parameters['tortuosity']
        vcl = parameters['viscous characteristic length']
        tcl = parameters['thermal characteristic length']
        air_gap = parameters['air gap']

        layer = dummy_struct.Add_JCA_Layer(t,fr,phi,tort,vcl,tcl)
        air = dummy_struct.Add_Air_Layer(thickness = air_gap)
        dummy_struct.frequency = self.frequency
        
        if self.opt_type == 'Gap' or self.opt_type == 'No Gap':
            s = dummy_struct.assemble_structure(layer,air)
            predicted = dummy_struct.absorption(s)[:,1]
            return(predicted)
        
        elif self.opt_type == 'Dual':
            no_gap_s = dummy_struct.assemble_structure(layer)
            gap_s = dummy_struct.assemble_structure(layer,air)
            
            no_gap_pred = dummy_struct.absorption(no_gap_s)[:,1]
            gap_pred = dummy_struct.absorption(gap_s)[:,1]
            
            return([no_gap_pred,gap_pred])
           
    def _error(self,
               x: list) -> float:
        """
        Function that is minimized in the optimization routine. Calculates the error between the measured (impedance tube)
        and predicted (TMM) absorption coefficients.

        Parameters
        ----------
        x (list):
            list structure containing the identified thickness, flow resistivity, porosity, tortuosity, viscous characteristic length,
            thermal characteristic length, and air gap thickness of the sample (in that order).

        Returns
        -------
        err (float):
            if opt_type is 'No Gap' or 'Gap' --> sum of the absolute square difference between measured and predicted absorption coefficients across all specified frequencies.
            if opt_type is 'Dual', the error for each mounting condition is averaged into a single error metric
        
        """

        t = x[0]
        fr = x[1]
        phi = x[2]
        tort = x[3]
        vcl = x[4]
        tcl = x[5]
        air = x[6]

        params = {'thickness':t,'flow resistivity':fr,'porosity':phi,'tortuosity':tort,'viscous characteristic length':vcl,'thermal characteristic length':tcl,'air gap':air}

        A = self._predictionJCA(params)
        
        if self.opt_type == 'No Gap':
            err = np.sum(np.abs(np.diff(A-self.meas_abs[0]))**2)
            return(err)
        
        elif self.opt_type == 'Gap':
            err = np.sum(np.abs(np.diff(A-self.meas_abs[1]))**2)
            return(err)
        
        elif self.opt_type == 'Dual':
            err = (np.sum(np.abs(np.diff(A[0]-self.meas_abs[0]))**2)+np.sum(np.abs(np.diff(A[1]-self.meas_abs[1]))**2))/2
            return(err)

    def _bounds(self,
                tort: float) -> tuple:
        """
        Defines the lower and upper boundary values for the parameters in the optimization routine. Bounds for thickness, flow resistivity,
        porosity, and air gap are calculated using the supplied "known" values and the uncertainty.  Tortuosity, viscous, and thermal characteristic
        length bounds are defined in:

        Atalla, Youssef & Panneton, R.. (2005). Inverse acoustical characterization of open cell porous 
        media using impedance tube measurements. Canadian Acoustics - Acoustique Canadienne. 33.

        Parameters
        ----------
        tort (float):
            tortuosity --> defined either as the initial guess or the most recent value returned by the optimization routine.

        Returns
        -------
        bounds (tuple):
            contains tuples of the lower and upper bounds for each parameter.
        
        """

        l_unc = 1-self.uncertainty
        u_unc = 1+self.uncertainty
        
        t_lb = l_unc*self.thickness
        t_ub = u_unc*self.thickness
        
        
        fr_lb = l_unc*self.flow_resistivity
        fr_ub = u_unc*self.flow_resistivity
        
        phi_lb = l_unc*self.porosity
        phi_ub = u_unc*self.porosity
        
        air_lb = l_unc*self.air_gap
        air_ub = u_unc*self.air_gap
        
        tort_lb = 1
        tort_ub = 4

        vcl_lb = ((1/3.3)*np.sqrt((8*tort*self.viscosity_temp)/(self.flow_resistivity*self.porosity)))/(1e-6)
        vcl_ub =((1/.3)*np.sqrt((8*tort*self.viscosity_temp)/(self.flow_resistivity*self.porosity)))/(1e-6)

        tcl_lb = vcl_lb
        tcl_ub = vcl_ub

        tort_b = (tort_lb,tort_ub)
        vcl_b = (vcl_lb,vcl_ub)
        tcl_b = (tcl_lb,tcl_ub)
        phi_b = (phi_lb,phi_ub)
        fr_b = (fr_lb,fr_ub)
        t_b = (t_lb,t_ub)
        air_b = (air_lb,air_ub)
        
        bounds = (t_b,fr_b,phi_b,tort_b,vcl_b,tcl_b,air_b)
        
        return (bounds)

    def Inverse(self,
                thickness: float,
                flow_resistivity: float,
                porosity: float,
                air_gap: float=0,
                uncertainty: float=0.01,
                early_stopping: float=1e-15,
                verbose: bool=False) -> dict:
        """
        Optimization routine for identifying the hard-to-measure JCA parameters (tortuosity, viscous, and thermal characteristic lengths).
        The routine uses the Sequential Least Squares Programming method (https://docs.scipy.org/doc/scipy/reference/optimize.minimize-slsqp.html)
        to minimize the error between predicted and actual absorption coefficients.

        If the early stopping criterion is not met after an initial guess, a unique grid search of the parameter space is crafted to help ensure
        the global minimum is found (ie the correct values for the parameters are identified).

        Parameters
        ----------
        thickness (float):
            The measured thickness of the sample [mm]
        
        flow_resistivity (float):
            The measured flow resistivity of the sample [Pa*s/m2]

        porosity (float):
            The measured porosity of the sample [-]

        air_gap (float):
            The measured impedance tube air gap behind the sample, if 'Gap' or 'Dual' mounting conditions are specified [mm]

        uncertainty (float):
            A measure of how uncertain the user is in the thickness, flow resistivity, porosity, and air gap measurements
            of the sample [0 - 1].  Increasing the uncertainty value will result in a wider search of the parameter space.

        early_stopping (float):
            criterion for ending the search early.  If the calculated error at any given step is less than this value, the routine will terminate
            and return the results. The default value of 1e-15 has been tested on a number of simulated cases.

        verbose (bool):
             If true, the progress of the optimization routine will print to the console.

        Returns
        -------
        result_dict (dict):
            dictionary containing the identified parameters associated with the lowest calculated error.
        
        """
        self.thickness = thickness
        self.flow_resistivity = flow_resistivity
        self.porosity = porosity
        self.air_gap = air_gap

        if uncertainty >= 0 and uncertainty <= 1:
            self.uncertainty = uncertainty

        else:
            print("Uncertainty must be between 0 and 1, reverting to default uncertainty of 1.0%!")
            self.uncertainty = 0.01

        init_tort = 2.5
        
        bnds = self._bounds(init_tort)
        x0 = [self.thickness,self.flow_resistivity,self.porosity,init_tort,(bnds[4][1]-bnds[4][0])/2,(bnds[5][1]-bnds[5][0])/2,self.air_gap] 
        cons = ({'type':'ineq','fun':lambda x:x[5]-x[4]})
        res = minimize(self._error,x0,method='SLSQP',bounds=bnds,constraints=cons,tol = 1e-50,options={'ftol':1e-50, 'maxiter':1000})
        err = res.fun
        results = np.around(res.x,decimals=3)

        if err < early_stopping:
            result_dict = {'thickness':results[0],'flow resistivity':results[1],'porosity':results[2],'tortuosity':results[3],'viscous characteristic length':results[4],'thermal characteristic length':results[5],'air gap':results[6],'error':err}  
            if verbose == True:
                print(f"Early stopping criteria has been met. The current lowest error is: {err}")
            return(result_dict)

        elif err > early_stopping:
            bnds_temp = self._bounds(results[3])
            tort_grid = [bnds_temp[3][0],(.75*bnds_temp[3][0]+.25*bnds_temp[3][1]),(.5*bnds_temp[3][0]+.5*bnds_temp[3][1]),(.25*bnds_temp[3][0]+.75*bnds_temp[3][1]),bnds_temp[3][1]]
            vcl_grid = [bnds_temp[4][0],(.75*bnds_temp[4][0]+.25*bnds_temp[4][1]),(.5*bnds_temp[4][0]+.5*bnds_temp[4][1]),(.25*bnds_temp[4][0]+.75*bnds_temp[4][1]),bnds_temp[4][1]]
            tcl_grid = [bnds_temp[5][0],(.75*bnds_temp[5][0]+.25*bnds_temp[5][1]),(.5*bnds_temp[5][0]+.5*bnds_temp[5][1]),(.25*bnds_temp[5][0]+.75*bnds_temp[5][1]),bnds_temp[5][1]]
            
            loop_len = len(tort_grid)*len(vcl_grid)*len(tcl_grid)
            i = 1
            
            for to in tort_grid: 
                for v in vcl_grid:
                    for tc in tcl_grid:
                        phi = self.porosity
                        fr = self.flow_resistivity
                        t = self.thickness
                        air = self.air_gap
                        x0 = [t,fr,phi,to,v,tc,air]

                        grid_search = round((i/loop_len)*100,2)
                        bnds = self._bounds(results[3])
                        res2 = minimize(self._error,x0,method='SLSQP',bounds=bnds,constraints=cons,tol = 1e-50,options={'ftol':1e-50, 'maxiter':1000})
                        
                        if verbose == True:
                            print(f"{grid_search}% of the parameter space has been searched. The current lowest error is: {err}")
                            

                        if res2.fun < err:
                            err = res2.fun
                            results = np.around(res2.x,decimals=3)
                           
                            if err < early_stopping:
                                result_dict = {'thickness':results[0],'flow resistivity':results[1],'porosity':results[2],'tortuosity':results[3],'viscous characteristic length':results[4],'thermal characteristic length':results[5],'air gap':results[6],'error':err}
                                if verbose == True:
                                    print(f"Early stopping criteria has been met. The current lowest error is: {err}")
                                return(result_dict)
                        i+=1
                       

            print("Warning: stopping criterion not met during parameter search. Double check known inputs and/or consider increasing the uncertainty value of the knowns.")
            result_dict = {'thickness': results[0],
                           'flow resistivity': results[1],
                           'porosity': results[2],
                           'tortuosity': results[3],
                           'viscous characteristic length': results[4],
                           'thermal characteristic length': results[5],
                           'air gap': results[6],
                           'error': err}
            
            return(result_dict)
        
    def Indirect(self,
                 thickness: float,
                 porosity: float,
                 flow_resistivity: float=None,
                 air_gap: float=0) -> dict:
        """
        Indirect method for identifying the hard-to-measure JCA parameters (tortuosity, viscous, and thermal characteristic lengths), based on:

        Panneton, R. & Salissou, Yacoubou. (2009). Indirect acoustical characterization of sound absorbing materials.. 
        The Journal of the Acoustical Society of America. 126. 2297. 10.1121/1.3249416.

        This method requires 'Dual' mounting conditions and either reflection or surface impedance measurements of the samples and 
        it is susceptible to uncertainty in the measurements of thickness, porosity, flow resistivity, air gap, 
        and the acoustic indicator(s) -- but the advantage is that it does not require a flow resistivity measurement in order to 
        estimate the JCA parameters.

        Parameters
        ----------
        thickness (float):
            The measured thickness of the sample [mm]

        porosity (float):
            The measured porosity of the sample [-]
        
        flow_resistivity (float):
            Optional, the measured flow resistivity of the sample [Pa*s/m2].  If no flow resistivity is specified, the value
            will be estimated automatically.

        air_gap (float):
            The measured impedance tube air gap behind the sample, if 'Gap' or 'Dual' mounting conditions are specified [mm]

        Returns
        -------
        result_dict (dict):
            dictionary containing the identified parameters associated with the lowest calculated error.
        
        """

        self.thickness = thickness
        self.porosity = porosity
        self.air_gap = air_gap
        
        if flow_resistivity == None:
            self.flow_resistivity = None
            return_preds = True
        else:
            self.flow_resistivity = flow_resistivity
            return_preds = False
        
        air_gap = self.air_gap/1000
        thickness = self.thickness/1000
        w = 2*np.pi*self.frequency
        k0 = w / self.soundspeed_temp
        
        if self.input_type == 'absorption':
            raise ValueError('Absorption data cannot be used for Indirect characterizations, please specify reflection or surface impedance data')
    
        elif self.opt_type == 'Gap' or self.opt_type == 'No Gap':
            raise ValueError('Dual mount types are required for Indirect characterizations')
        
        elif self.input_type == 'surface':
            Zs_NG = self.no_gap_data[:,1]
            Zs_G = self.gap_data[:,1]

        elif self.input_type == 'reflection':
            Zs_NG = self.Z0*((1+self.no_gap_data[:,1])/(1-self.no_gap_data[:,1]))
            Zs_G = self.Z0*((1+self.gap_data[:,1])/(1-self.gap_data[:,1]))
            

        T11A = np.cos(k0*air_gap)
        T21A = (1j/self.Z0)*np.sin(k0*air_gap)
        Zs_A = T11A/T21A

        Zp = np.sqrt((Zs_NG*(Zs_G-Zs_A))+(Zs_G*Zs_A))
        kp = np.arctan(Zp/(1j*Zs_NG))/thickness
        
        test = np.column_stack((self.frequency,np.real(Zp/(1j*Zs_NG))))
        
        try:
            cutoff = test[np.where(test[:-1] * test[1:] < 0 )[0]]
            cutoff = np.abs(np.min(cutoff[:,0]))
        except ValueError:
            cutoff = np.max(self.frequency)
        
        peff = Zp*kp/w
        keff = w*np.divide(Zp,kp)
        
        re_peff = np.real(peff)
        im_peff = np.imag(peff)
        
        fr = np.column_stack((self.frequency**2,-im_peff*w))
        fr_curve = fr[np.where(fr[:,0] <= cutoff**2)]
        slope1,intercept1,r_value1,p_value1,std_err1 = scipy.stats.linregress(fr_curve[:,0],fr_curve[:,1])
        
        jca_fr = np.abs(intercept1)
        
        #######################################################################################################
        
        if return_preds == False:
            tort = (self.porosity/self.density_temp)*(re_peff-np.sqrt((im_peff**2)-((self.flow_resistivity/w)**2)))
        elif return_preds == True:
            tort = (self.porosity/self.density_temp)*(re_peff-np.sqrt((im_peff**2)-((jca_fr/w)**2)))
        tort1 = np.column_stack((self.frequency,tort))
        tort_curve = tort1[np.where(tort1[:,0] <= cutoff )]
        jca_tort = np.abs(np.mean(tort_curve[:,1]))
        
        #######################################################################################################
        
        phi_num = (self.density_temp*jca_tort)
        if return_preds == False:
            phi_denom = re_peff-np.sqrt((im_peff**2)-((self.flow_resistivity/w)**2))
        elif return_preds == True:
            phi_denom = re_peff-np.sqrt((im_peff**2)-((jca_fr/w)**2))
        phi = phi_num/phi_denom
        phi1 = np.column_stack((self.frequency,phi))
        phi_curve = phi1[np.where((phi1[:,0] <= cutoff) & (phi1[:,0] >= 500))]
        jca_phi = np.abs(np.mean(phi_curve[:,1]))

        ####################################################################################################### 
        
        if return_preds == False:
            vcl = jca_tort*np.sqrt((2*self.density_temp*self.viscosity_temp)/((w*self.porosity*im_peff)*((self.density_temp*jca_tort)-(self.porosity*re_peff))))/(1e-6)
        elif return_preds == True:
            vcl = jca_tort*np.sqrt((2*self.density_temp*self.viscosity_temp)/((w*jca_phi*im_peff)*((self.density_temp*jca_tort)-(jca_phi*re_peff))))/(1e-6)
        vcl1 = np.column_stack((self.frequency,vcl))
        vcl_curve = vcl1[np.where(vcl1[:,0] <= cutoff )]
        jca_vcl = np.abs(np.mean(vcl_curve[:,1]))
        
        #######################################################################################################        
        
        tcl0 = np.sqrt((2*self.viscosity_temp)/(w*self.density_temp*self.Pr_temp))
        if return_preds == False:
            tcl1 = ((1-((keff*self.porosity)/(self.gamma_temp*self.P0)))/(1-((keff*self.porosity)/self.P0)))**2
        elif return_preds == True:
            tcl1 = ((1-((keff*jca_phi)/(self.gamma_temp*self.P0)))/(1-((keff*jca_phi)/self.P0)))**2
        tcl2 = 1/(-np.imag(tcl1))
        tcl3 = np.sqrt(2*tcl2)       
        tcl4 = tcl0*tcl3/(1e-6)
        tcl5 = np.column_stack((self.frequency,tcl4))
        tcl_curve = tcl5[np.where((tcl5[:,0] <= cutoff) & (tcl5[:,0] >= 500))]
        jca_tcl = np.abs(np.mean(tcl_curve[:,1]))
        
        #######################################################################################################           
        
        if return_preds == False:
            kn0 = (self.porosity*self.viscosity_temp)/(w*self.density_temp*self.Pr_temp)
        elif return_preds == True:
            kn0 = (jca_phi*self.viscosity_temp)/(w*self.density_temp*self.Pr_temp)
        kn1 = 1/np.sqrt(-np.real(tcl1))
        kn3 = kn0*kn1/(1e-10)
        kn4 = np.column_stack((self.frequency,kn3))
        k0_curve = kn4[np.where((kn4[:,0] <= cutoff) & (kn4[:,0] >= 500))]
        jca_k0 = np.abs(np.mean(k0_curve[:,1]))
        
        result_dict = {'thickness': self.thickness,
                       'flow resistivity': self.flow_resistivity,
                       'porosity': self.porosity,
                       'tortuosity': round(jca_tort,4),
                        'viscous characteristic length': round(jca_vcl,4),
                        'thermal characteristic length': round(jca_tcl,4),
                        'air gap': self.air_gap} 
        
        full_predicted_dict = {'thickness': self.thickness,
                               'flow resistivity': jca_fr,
                               'porosity': jca_phi,
                               'tortuosity': jca_tort,
                               'viscous characteristic length': jca_vcl,
                               'thermal characteristic length': jca_tcl,
                               'air gap': self.air_gap} 
        
        if return_preds == False:
            return(result_dict)

        elif return_preds == True:
            return(full_predicted_dict)

    
    def Hybrid(self,
               thickness: float,
               porosity: float,
               flow_resistivity: float=None,
               air_gap: float=0,
               uncertainty: float=0.01,
               early_stopping: float=1e-15,
               verbose: bool=False) -> dict:
        """
        The Hybrid routine uses both the inverse and indirect characterization methods to identify the JCA parameters.  The parameters are 
        first estimated using the indirect method and the error between measured and estimated absorption coefficients is determined.

        If the early stopping criterion is not met using the indirect method, the estimate is then used as the initial guess for the inverse
        procedure and to calculate the bounds of the grid search.
        
        This method requires 'Dual' mounting conditions and either reflection or surface impedance measurements of the samples.
        The advantage of using this procedure compared to the inverse or indirect methods alone are:
            Inverse: The Hybrid method does not require flow resitivity to be known.
            Indirect:  The Hybrid method is much less susceptible to uncertainty in the measurements.

        Parameters
        ----------
        thickness (float):
            The measured thickness of the sample [mm]
        
        flow_resistivity (float):
            The measured flow resistivity of the sample [Pa*s/m2]

        porosity (float):
            The measured porosity of the sample [-]

        air_gap (float):
            The measured impedance tube air gap behind the sample, if 'Gap' or 'Dual' mounting conditions are specified [mm]

        uncertainty (float):
            A measure of how uncertain the user is in the thickness, flow resistivity, porosity, and air gap measurements
            of the sample [0 - 1].  Increasing the uncertainty value will result in a wider search of the parameter space.

        early_stopping (float):
            criterion for ending the search early.  If the calculated error at any given step is less than this value, the routine will terminate
            and return the results. The default value of 1e-15 has been tested on a number of simulated cases.

        verbose (bool):
             If true, the progress of the optimization routine will print to the console.

        Returns
        -------
        result_dict (dict):
            dictionary containing the identified parameters associated with the lowest calculated error.
        
        """
        self.thickness = thickness
        self.porosity = porosity
        self.air_gap = air_gap
        if uncertainty >= 0 and uncertainty <= 1:
            self.uncertainty = uncertainty

        else:
            print("Uncertainty must be between 0 and 1, reverting to default uncertainty of 1.0%!")
            self.uncertainty = 0.01
        

        indirect_results = self.Indirect(self.thickness,self.porosity,flow_resistivity,self.air_gap)
        
        if np.isnan(indirect_results['thermal characteristic length']):
            indirect_results['thermal characteristic length'] = indirect_results['viscous characteristic length']
        
        self.flow_resistivity = indirect_results['flow resistivity']

        init_tort = indirect_results['tortuosity']
        
        bnds = self._bounds(init_tort)
        x0 = [self.thickness,indirect_results['flow resistivity'],indirect_results['porosity'],init_tort,indirect_results['viscous characteristic length'],indirect_results['thermal characteristic length'],self.air_gap] 
        cons = ({'type':'ineq','fun':lambda x:x[5]-x[4]})
        res = minimize(self._error,x0,method='SLSQP',bounds=bnds,constraints=cons,tol = 1e-50,options={'ftol':1e-50, 'maxiter':1000})
        err = res.fun
        results = np.around(res.x,decimals=3)

        if err < early_stopping:
            result_dict = {'thickness':results[0],'flow resistivity':results[1],'porosity':results[2],'tortuosity':results[3],'viscous characteristic length':results[4],'thermal characteristic length':results[5],'air gap':results[6],'error':err}  
            if verbose == True:
                print(f"Early stopping criteria has been met. The current lowest error is: {err}")
            return(result_dict)

        elif err > early_stopping:
            bnds_temp = self._bounds(results[3])
            tort_grid = [bnds_temp[3][0],(.75*bnds_temp[3][0]+.25*bnds_temp[3][1]),(.5*bnds_temp[3][0]+.5*bnds_temp[3][1]),(.25*bnds_temp[3][0]+.75*bnds_temp[3][1]),bnds_temp[3][1]]
            vcl_grid = [bnds_temp[4][0],(.75*bnds_temp[4][0]+.25*bnds_temp[4][1]),(.5*bnds_temp[4][0]+.5*bnds_temp[4][1]),(.25*bnds_temp[4][0]+.75*bnds_temp[4][1]),bnds_temp[4][1]]
            tcl_grid = [bnds_temp[5][0],(.75*bnds_temp[5][0]+.25*bnds_temp[5][1]),(.5*bnds_temp[5][0]+.5*bnds_temp[5][1]),(.25*bnds_temp[5][0]+.75*bnds_temp[5][1]),bnds_temp[5][1]]
            
            loop_len = len(tort_grid)*len(vcl_grid)*len(tcl_grid)
            i = 1
            
            for to in tort_grid: 
                for v in vcl_grid:
                    for tc in tcl_grid:
                        phi = self.porosity
                        fr = self.flow_resistivity
                        t = self.thickness
                        air = self.air_gap
                        x0 = [t,fr,phi,to,v,tc,air]

                        grid_search = round((i/loop_len)*100,2)
                        bnds = self._bounds(results[3])
                        res2 = minimize(self._error,x0,method='SLSQP',bounds=bnds,constraints=cons,tol = 1e-50,options={'ftol':1e-50, 'maxiter':1000})
                        
                        if verbose == True:
                            print(f"{grid_search}% of the parameter space has been searched. The current lowest error is: {err}")
                            

                        if res2.fun < err:
                            err = res2.fun
                            results = np.around(res2.x,decimals=3)
                           
                            if err < early_stopping:
                                result_dict = {'thickness':results[0],'flow resistivity':results[1],'porosity':results[2],'tortuosity':results[3],'viscous characteristic length':results[4],'thermal characteristic length':results[5],'air gap':results[6],'error':err}
                                if verbose == True:
                                    print(f"Early stopping criteria has been met. The current lowest error is: {err}")
                                return(result_dict)
                        i+=1
                       

            print("Warning: stopping criterion not met during parameter search. Double check known inputs and/or consider increasing the uncertainty value of the knowns.")
            result_dict = {'thickness': results[0],
                           'flow resistivity': results[1],
                           'porosity': results[2],
                           'tortuosity': results[3],
                           'viscous characteristic length': results[4],
                           'thermal characteristic length': results[5],
                           'air gap': results[6],
                           'error': err}
            
            return(result_dict)
        
        
    def stats(self,
              parameters:dict) -> dict:
        """
        Calculates statistics about the parameters identified in the optimization routine via linear regression of the predicted vs measured
        absorption coefficients.

        Parameters
        ----------
        parameters (dict):
            dictionary containing the identified thickness, flow resistivity, porosity, tortuosity, 
            viscous characteristic length, thermal characteristic length, and air gap of the sample.

        Returns
        -------
        stats (dict):
            dictionary containing the slope, intercept, r value, p value, and std error returned from the linear regression.
            If 'Dual' opt_type is specified, the statistics for each mounting condition are averaged.
        
        """
        if self.opt_type == 'No Gap':
            x = self.meas_abs[0]
            y = self._predictionJCA(parameters)
            
            slope,intercept,r_value,p_value,std_err = scipy.stats.linregress(x,y)
            stats = {'slope':slope,'intercept':intercept,'r_value':r_value,'p_value':p_value,'std_err':std_err}
            
            return(stats)
        
        elif self.opt_type == 'Gap':
            x = self.meas_abs[1]
            y = self._predictionJCA(parameters)
            
            slope,intercept,r_value,p_value,std_err = scipy.stats.linregress(x,y)
            stats = {'slope':slope,'intercept':intercept,'r_value':r_value,'p_value':p_value,'std_err':std_err}

            return(stats)
        
        elif self.opt_type == 'Dual':
            x1 = self.meas_abs[0]
            y1 = self._predictionJCA(parameters)[0]
            
            slope1,intercept1,r_value1,p_value1,std_err1 = scipy.stats.linregress(x1,y1)
            
            x2 = self.meas_abs[1]
            y2 = self._predictionJCA(parameters)[1]
            
            slope2,intercept2,r_value2,p_value2,std_err2 = scipy.stats.linregress(x2,y2)
            
            slope = np.abs((slope1+slope2)/2)
            intercept = np.abs((intercept1+intercept2)/2)
            r_value = np.abs((r_value1+r_value2)/2)
            p_value = np.abs((p_value1+p_value2)/2)
            std_err = np.abs((std_err1+std_err2)/2)
                        
            stats = {'slope':slope,'intercept':intercept,'r_value':r_value,'p_value':p_value,'std_err':std_err}
            
            return(stats)

    def plot_comparison(self,
                        parameters:dict) -> None:
        """
        Plots the predicted and measured frequency dependent absorption coefficients.

        Parameters
        ----------
        parameters (dict):
            dictionary containing the identified thickness, flow resistivity, porosity, tortuosity, 
            viscous characteristic length, thermal characteristic length, and air gap of the sample.
        
        """

        f, ax = plt.subplots(1)
        
        if self.opt_type == 'No Gap':
            actual = self.meas_abs[0]
            predicted = self._predictionJCA(parameters)
            
            ax.plot(self.frequency,actual,label='Actual')
            ax.plot(self.frequency,predicted,label='Predicted')

            ax.legend(loc="lower right")
            ax.set_ylim(bottom=0)
    
            plt.show()

        
        elif self.opt_type == 'Gap':
            actual = self.meas_abs[1]
            predicted = self._predictionJCA(parameters)
            
            ax.plot(self.frequency,actual,label='Actual')
            ax.plot(self.frequency,predicted,label='Predicted')

            ax.legend(loc="lower right")
            ax.set_ylim(bottom=0)
    
            plt.show()
        
        elif self.opt_type == 'Dual':
            no_gap_pred = self._predictionJCA(parameters)[0]
            gap_pred = self._predictionJCA(parameters)[1]
            no_gap_actual = self.meas_abs[0]
            gap_actual = self.meas_abs[1]
    
            ax.plot(self.frequency,no_gap_actual,label='No Gap Actual')
            ax.plot(self.frequency,no_gap_pred,label='No Gap Predicted')
            ax.plot(self.frequency,gap_actual,label='Gap Actual')
            ax.plot(self.frequency,gap_pred,label='Gap Predicted')
            
            ax.legend(loc="lower right")
            ax.set_ylim(bottom=0)
    
            plt.show()

    def to_csv(self,
               FileName: str,
               parameters: dict) -> None:
        """
        Saves the identified parameters and the measured/predicted absorption curves to a csv file.

        Parameters
        ----------
        FileName (str):
            Name of the csv file to save data to.

        parameters (dict):
            dictionary containing the identified thickness, flow resistivity, porosity, tortuosity, 
            viscous characteristic length, thermal characteristic length, and air gap of the sample.
        
        """
        
        if ".csv" in FileName:
            file = FileName
        else:
            file = FileName+".csv"


        if self.opt_type == 'No Gap':
            actual = self.meas_abs[0]
            predicted = self._predictionJCA(parameters)
            all_data = {'frequency':self.frequency,'measured':actual,'predicted':predicted}

        
        elif self.opt_type == 'Gap':
            actual = self.meas_abs[1]
            predicted = self._predictionJCA(parameters)
            all_data = {'frequency':self.frequency,'measured':actual,'predicted':predicted}
        
        elif self.opt_type == 'Dual':
            no_gap_pred = self._predictionJCA(parameters)[0]
            gap_pred = self._predictionJCA(parameters)[1]
            no_gap_actual = self.meas_abs[0]
            gap_actual = self.meas_abs[1]
        
            all_data = {'frequency':self.frequency,'no gap measured':no_gap_actual,'no gap predicted':no_gap_pred,'gap measured':gap_actual,'gap predicted':gap_pred}

        with open(file, 'w') as f:
            [f.write('{0},{1}\n'.format(key, value)) for key, value in parameters.items()]
            writer = csv.writer(f, delimiter = ",",lineterminator='\n')
            writer.writerow(all_data.keys())
            writer.writerows(zip(*all_data.values()))

    def load_to_array(self,
                      FileName: str,
                      type: str ='complex') -> None:
        """
        Loads data from csv or excel file.

        Parameters
        ----------
        FileName (str):
            Name of the file to load data from

        type (str):
            type of data being loaded -- either complex or floating point data
        
        """
        if type == 'complex':
            try:
                data = np.asarray(pd.read_csv(FileName,header=None).map(lambda s: np.complex128(s.replace('i', 'j'))))
            except Exception:
                data = np.asarray(pd.read_excel(FileName,header=None).map(lambda s: np.complex128(s.replace('i', 'j'))))

        elif type == 'float':
            try:
                data = np.asarray(pd.read_csv(FileName,header=None))
            except Exception:
                data = np.asarray(pd.read_excel(FileName,header=None))  
        
        return (data)
            
    def to_database(self,
                    parameters: dict,
                    layer_name: str) -> None:
        """
        Saves the identified parameters as a new layer in a database.

        Parameters
        ----------
        parameters (dict):
            dictionary containing the identified thickness, flow resistivity, porosity, tortuosity, 
            viscous characteristic length, thermal characteristic length, and air gap of the sample.

        layer_name (str):
            Specifies the name of the layer.  Must be a unique identifier.
        """

        thickness = parameters['thickness']
        flow_resistivity = parameters['flow resistivity']
        porosity = parameters['porosity']
        tortuosity = parameters['tortuosity']
        viscous_characteristic_length = parameters['viscous characteristic length']
        thermal_characteristic_length = parameters['thermal characteristic length']
        
        s = AcoustiBase()
        data = s.pull('LAYER')
        id1 = len(data)+1
        params = [id1,layer_name,'identified','null',thickness,flow_resistivity,porosity,tortuosity,viscous_characteristic_length,thermal_characteristic_length,'null','null','null','null','null','null']
        s.execute(params,'LAYER')
        s.commit()
        s.close()