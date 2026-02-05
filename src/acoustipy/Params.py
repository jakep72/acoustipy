import numpy as np
import pandas as pd
import torch
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
    # Valid options for mount_type and input_type
    VALID_MOUNT_TYPES = ('No Gap', 'Gap', 'Dual')
    VALID_INPUT_TYPES = ('absorption', 'reflection', 'surface')

    def __init__(self,
                 mount_type: str = "No Gap",
                 no_gap_file: str = None,
                 gap_file: str = None,
                 input_type: str = 'absorption',
                 air_temperature: float = None,
                 sound_speed: float = 343.152,
                 air_density: float = 1.2058,
                 Cp: float = 1.004425,
                 Cv: float = 0.717425,
                 viscosity: float = 1.825e-05,
                 Pr: float = 0.7157,
                 P0: float = 101325,
                 device: str = 'cpu'
                 ):
        """
        Initialize AcousticID parameter identification object.
        
        Parameters
        ----------
        device : str, optional
            Device for computations ('cpu' or 'cuda'). Default is 'cpu'.
            GPU acceleration is primarily used by the ML method.
        """
        # Validate mount_type
        if mount_type not in self.VALID_MOUNT_TYPES:
            raise ValueError(
                f"Invalid mount_type '{mount_type}'. "
                f"Must be one of: {', '.join(self.VALID_MOUNT_TYPES)}"
            )
        
        # Validate input_type
        if input_type not in self.VALID_INPUT_TYPES:
            raise ValueError(
                f"Invalid input_type '{input_type}'. "
                f"Must be one of: {', '.join(self.VALID_INPUT_TYPES)}"
            )
        
        # Validate required files based on mount_type
        self._validate_files(mount_type, no_gap_file, gap_file)

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
        self.device = device
        self.thickness = None
        self.flow_resistivity = None
        self.porosity = None
        self.uncertainty = None
        self.air_gap = None

        # Load data based on mount type and input type
        data_type = 'float' if input_type == 'absorption' else 'complex'
        
        if mount_type == 'No Gap':
            self.no_gap_data = self.load_to_array(no_gap_file, type=data_type)

        elif mount_type == 'Gap':
            self.gap_data = self.load_to_array(gap_file, type=data_type)

        elif mount_type == 'Dual':
            self.no_gap_data = self.load_to_array(no_gap_file, type=data_type)
            self.gap_data = self.load_to_array(gap_file, type=data_type)
    
    @staticmethod
    def _validate_files(mount_type: str, no_gap_file: str, gap_file: str) -> None:
        """
        Validate that required files are provided and exist.
        
        Parameters
        ----------
        mount_type : str
            The mounting type being used.
        no_gap_file : str
            Path to the no-gap data file.
        gap_file : str
            Path to the gap data file.
            
        Raises
        ------
        ValueError
            If required files are not provided or do not exist.
        """
        import os
        
        if mount_type in ('No Gap', 'Dual'):
            if no_gap_file is None:
                raise ValueError(f"no_gap_file is required for mount_type='{mount_type}'")
            if not os.path.exists(no_gap_file):
                raise FileNotFoundError(f"no_gap_file not found: {no_gap_file}")
        
        if mount_type in ('Gap', 'Dual'):
            if gap_file is None:
                raise ValueError(f"gap_file is required for mount_type='{mount_type}'")
            if not os.path.exists(gap_file):
                raise FileNotFoundError(f"gap_file not found: {gap_file}")
    
    @staticmethod
    def _validate_physical_params(thickness: float, flow_resistivity: float, 
                                   porosity: float, air_gap: float = 0) -> None:
        """
        Validate that material parameters are physically reasonable.
        
        Parameters
        ----------
        thickness : float
            Sample thickness in mm.
        flow_resistivity : float
            Flow resistivity in Pa*s/m^2.
        porosity : float
            Open porosity (0-1).
        air_gap : float, optional
            Air gap thickness in mm.
            
        Raises
        ------
        ValueError
            If parameters are physically invalid.
        """
        if thickness <= 0:
            raise ValueError(f"Thickness must be positive, got {thickness}")
        
        if flow_resistivity <= 0:
            raise ValueError(f"Flow resistivity must be positive, got {flow_resistivity}")
        
        if not (0 < porosity <= 1):
            raise ValueError(f"Porosity must be in range (0, 1], got {porosity}")
        
        if air_gap < 0:
            raise ValueError(f"Air gap cannot be negative, got {air_gap}")
        
    @property
    def frequency(self):
        #frequency range of interest
        
        if self.opt_type == 'No Gap':
            no_gap_freq = torch.tensor(self.no_gap_data[:,0])
            return(no_gap_freq)
        
        elif self.opt_type == 'Gap':
            gap_freq = self.gap_data[:,0]
            return(gap_freq)
        
        elif self.opt_type == 'Dual':
            no_gap_freq = torch.tensor(self.no_gap_data[:, 0]).real
            gap_freq = torch.tensor(self.gap_data[:, 0]).real
        
            if not torch.equal(no_gap_freq, gap_freq):
                raise ValueError("Frequencies must match between no gap and gap absorption curves")
            else:
                return no_gap_freq


    @property
    def meas_abs(self):
        #measured absorption coefficients
        if self.input_type == 'absorption':
            if self.opt_type == 'No Gap':
                no_gap_abs = torch.tensor(self.no_gap_data[:,1])
                return([no_gap_abs,None])
            
            elif self.opt_type == 'Gap':
                gap_abs = torch.tensor(self.gap_data[:,1])
                return([None,gap_abs])
            
            elif self.opt_type == 'Dual':
                no_gap_abs = torch.tensor(self.no_gap_data[:,1])
                gap_abs = torch.tensor(self.gap_data[:,1])
                return([no_gap_abs,gap_abs])

        elif self.input_type == 'reflection':
            if self.opt_type == 'No Gap':
                no_gap_abs = torch.tensor(1-abs(self.no_gap_data[:,1])**2)
                return([no_gap_abs,None])
            
            elif self.opt_type == 'Gap':
                gap_abs = torch.tensor(1-abs(self.gap_data[:,1])**2)
                return([None,gap_abs])
            
            elif self.opt_type == 'Dual':
                no_gap_abs = torch.tensor(1-abs(self.no_gap_data[:,1])**2,dtype=torch.float32).real
                gap_abs = torch.tensor(1-abs(self.gap_data[:,1])**2,dtype=torch.float32).real
                return([no_gap_abs,gap_abs])
        
        elif self.input_type == 'surface':
            if self.opt_type == 'No Gap':
                ngzs = self.no_gap_data
                ngr = (ngzs-self.Z0)/(ngzs+self.Z0)
                no_gap_abs = torch.tensor(1-abs(ngr[:,1])**2)
                return([no_gap_abs,None])
            
            elif self.opt_type == 'Gap':
                gzs = self.gap_data
                gr = (gzs-self.Z0)/(gzs+self.Z0)
                gap_abs = torch.tensor(1-abs(gr[:,1])**2)
                return([None,gap_abs])
            
            elif self.opt_type == 'Dual':
                ngzs = self.no_gap_data
                ngr = (ngzs-self.Z0)/(ngzs+self.Z0)
                no_gap_abs = torch.tensor(1-abs(ngr[:,1])**2)

                gzs = self.gap_data
                gr = (gzs-self.Z0)/(gzs+self.Z0)
                gap_abs = torch.tensor(1-abs(gr[:,1])**2)
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
                       parameters: dict,
                       optimize: bool) -> np.ndarray:
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
        if optimize:
            t = parameters[0]
            fr = parameters[1]
            phi = parameters[2]
            tort = parameters[3]
            vcl = parameters[4]
            tcl = parameters[5]
            air_gap = parameters[6]
        else:
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
               x: list) -> tuple[float, np.ndarray]:
        """
        Function that is minimized in the optimization routine. Calculates the error between the measured (impedance tube)
        and predicted (TMM) absorption coefficients.

        Parameters
        ----------
        x : list
            List containing the identified thickness, flow resistivity, porosity, tortuosity, viscous characteristic length,
            thermal characteristic length, and air gap thickness of the sample (in that order).

        Returns
        -------
        err : float
            If opt_type is 'No Gap' or 'Gap': sum of the absolute square difference between measured and predicted 
            absorption coefficients across all specified frequencies.
            If opt_type is 'Dual': the error for each mounting condition is averaged into a single error metric.
        grad : np.ndarray
            Gradient of the error with respect to the parameters.
        """
        params = torch.tensor(x, requires_grad=True)

        A = self._predictionJCA(params, optimize=True)
        
        if self.opt_type == 'No Gap':
            err = np.sum(np.abs(np.diff(A-self.meas_abs[0]))**2)
        
        elif self.opt_type == 'Gap':
            err = np.sum(np.abs(np.diff(A-self.meas_abs[1]))**2)
        
        elif self.opt_type == 'Dual':
            gap_loss = ((torch.sub(A[0],self.meas_abs[0])**2)).sum()
            no_gap_loss = ((torch.sub(A[1],self.meas_abs[1])**2)).sum()
            err = torch.mean(gap_loss+no_gap_loss,dtype=torch.double)
        
        err.backward()

        return err.data.cpu().numpy(), params.grad.data.cpu().numpy()

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
                early_stopping: float=1e-10,
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
        result_dict : dict
            Dictionary containing the identified parameters associated with the lowest calculated error.
        
        Raises
        ------
        ValueError
            If input parameters are physically invalid.
        """
        # Validate physical constraints
        self._validate_physical_params(thickness, flow_resistivity, porosity, air_gap)
        
        self.thickness = thickness
        self.flow_resistivity = flow_resistivity
        self.porosity = porosity
        self.air_gap = air_gap

        if 0 <= uncertainty <= 1:
            self.uncertainty = uncertainty
        else:
            warnings.warn("Uncertainty must be between 0 and 1, reverting to default uncertainty of 1.0%!")
            self.uncertainty = 0.01

        init_tort = 2.5
        
        bnds = self._bounds(init_tort)
        
        x0 = torch.tensor([self.thickness,self.flow_resistivity,self.porosity,init_tort,(bnds[4][1]-bnds[4][0])/2,(bnds[5][1]-bnds[5][0])/2,self.air_gap])
        cons = ({'type':'ineq','fun':lambda x:x[5]-x[4]})
        res = minimize(self._error,x0,method='SLSQP',bounds=bnds,constraints=cons,jac=True, options={'ftol':1e-50, 'maxiter':5000})
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
                        res2 = minimize(self._error,x0,method='SLSQP',bounds=bnds,constraints=cons,jac=True,tol = 1e-50,options={'ftol':1e-50, 'maxiter':1000})
                        
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
        w = 2*np.pi*self.frequency.numpy()
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
        
        test = np.column_stack((self.frequency.numpy(),np.real(Zp/(1j*Zs_NG))))
        
        try:
            cutoff = test[np.where(test[:-1] * test[1:] < 0 )[0]]
            cutoff = np.abs(np.min(cutoff[:,0]))
        except ValueError:
            cutoff = np.max(self.frequency.numpy())
        
        peff = Zp*kp/w
        keff = w*np.divide(Zp,kp)
        
        re_peff = np.real(peff)
        im_peff = np.imag(peff)
        
        fr = np.column_stack((self.frequency.numpy()**2,-im_peff*w))
        fr_curve = fr[np.where(fr[:,0] <= cutoff**2)]
        slope1,intercept1,r_value1,p_value1,std_err1 = scipy.stats.linregress(fr_curve[:,0],fr_curve[:,1])
        
        jca_fr = np.abs(intercept1)
        
        
        #######################################################################################################
        
        if return_preds == False:
            tort = (self.porosity/self.density_temp)*(re_peff-np.sqrt((im_peff**2)-((self.flow_resistivity/w)**2)))
        elif return_preds == True:
            tort = (self.porosity/self.density_temp)*(re_peff-np.sqrt((im_peff**2)-((jca_fr/w)**2)))
        
        tort1 = np.column_stack((self.frequency.numpy(),tort))
        
        tort_curve = tort1[np.where(tort1[:,0] <= cutoff )][:,1]
        tort_curve_dropped = tort_curve[np.isfinite(tort_curve)]
        jca_tort = np.abs(np.mean(tort_curve_dropped))
        
        
        #######################################################################################################
        
        phi_num = (self.density_temp*jca_tort)
        if return_preds == False:
            phi_denom = re_peff-np.sqrt((im_peff**2)-((self.flow_resistivity/w)**2))
        elif return_preds == True:
            phi_denom = re_peff-np.sqrt((im_peff**2)-((jca_fr/w)**2))
        phi = phi_num/phi_denom
        phi1 = np.column_stack((self.frequency.numpy(),phi))
        phi_curve = phi1[np.where((phi1[:,0] <= cutoff) & (phi1[:,0] >= 500))]
        jca_phi = np.abs(np.mean(phi_curve[:,1]))

        ####################################################################################################### 
        
        if return_preds == False:
            vcl = jca_tort*np.sqrt((2*self.density_temp*self.viscosity_temp)/((w*self.porosity*im_peff)*((self.density_temp*jca_tort)-(self.porosity*re_peff))))/(1e-6)
        elif return_preds == True:
            vcl = jca_tort*np.sqrt((2*self.density_temp*self.viscosity_temp)/((w*jca_phi*im_peff)*((self.density_temp*jca_tort)-(jca_phi*re_peff))))/(1e-6)
        vcl1 = np.column_stack((self.frequency.numpy(),vcl))
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
        tcl5 = np.column_stack((self.frequency.numpy(),tcl4))
        tcl_curve = tcl5[np.where((tcl5[:,0] <= cutoff) & (tcl5[:,0] >= 500))]
        jca_tcl = np.abs(np.mean(tcl_curve[:,1]))
        
        #######################################################################################################           
        
        if return_preds == False:
            kn0 = (self.porosity*self.viscosity_temp)/(w*self.density_temp*self.Pr_temp)
        elif return_preds == True:
            kn0 = (jca_phi*self.viscosity_temp)/(w*self.density_temp*self.Pr_temp)
        kn1 = 1/np.sqrt(-np.real(tcl1))
        kn3 = kn0*kn1/(1e-10)
        kn4 = np.column_stack((self.frequency.numpy(),kn3))
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
               early_stopping: float=1e-10,
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
        result_dict : dict
            Dictionary containing the identified parameters associated with the lowest calculated error.
        
        Raises
        ------
        ValueError
            If input parameters are physically invalid.
        """
        # Validate physical constraints (flow_resistivity can be None for Hybrid)
        if thickness <= 0:
            raise ValueError(f"Thickness must be positive, got {thickness}")
        if not (0 < porosity <= 1):
            raise ValueError(f"Porosity must be in range (0, 1], got {porosity}")
        if air_gap < 0:
            raise ValueError(f"Air gap cannot be negative, got {air_gap}")
        if flow_resistivity is not None and flow_resistivity <= 0:
            raise ValueError(f"Flow resistivity must be positive, got {flow_resistivity}")
        
        self.thickness = thickness
        self.porosity = porosity
        self.air_gap = air_gap
        
        if 0 <= uncertainty <= 1:
            self.uncertainty = uncertainty
        else:
            warnings.warn("Uncertainty must be between 0 and 1, reverting to default uncertainty of 1.0%!")
            self.uncertainty = 0.01

        indirect_results = self.Indirect(thickness=self.thickness, porosity=self.porosity, 
                                         flow_resistivity=flow_resistivity, air_gap=self.air_gap)
        
        if np.isnan(indirect_results['thermal characteristic length']):
            indirect_results['thermal characteristic length'] = indirect_results['viscous characteristic length']
        
        self.flow_resistivity = indirect_results['flow resistivity']

        init_tort = indirect_results['tortuosity']
        
        bnds = self._bounds(init_tort)
        x0 = torch.tensor([self.thickness,indirect_results['flow resistivity'],indirect_results['porosity'],init_tort,indirect_results['viscous characteristic length'],indirect_results['thermal characteristic length'],self.air_gap])
        cons = ({'type':'ineq','fun':lambda x:x[5]-x[4]})
        res = minimize(self._error,x0,method='SLSQP',bounds=bnds,constraints=cons, jac=True, tol = 1e-50,options={'ftol':1e-50, 'maxiter':1000})
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
                        res2 = minimize(self._error,x0,method='SLSQP',bounds=bnds,constraints=cons,jac=True, tol = 1e-50,options={'ftol':1e-50, 'maxiter':1000})
                        
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
        
    def _criterion(self, y1: torch.Tensor, y2: torch.Tensor, params: dict) -> torch.Tensor:
        """
        Loss function for the ML optimization method with penalty terms for invalid parameters.

        Parameters
        ----------
        y1 : torch.Tensor
            Predicted absorption coefficients.
        y2 : torch.Tensor
            Measured absorption coefficients.
        params : dict
            Dictionary of current parameter values (normalized).

        Returns
        -------
        err : torch.Tensor
            Weighted error value including penalties for physically invalid parameters.
        """
        err = 1e12 * torch.sum(torch.diff(y1 - y2) ** 2)
        if params['vcl'] > params['tcl']:
            err = 2 * err
        if params['fr'] > 1 or params['fr'] < 0:
            err = 2 * err
        if params['phi'] > 1 or params['phi'] < 0.001:
            err = 10 * err
        if params['tau'] > 1 or params['tau'] < 0.2:
            err = 2 * err
        if params['vcl'] > 1 or params['tcl'] > 1:
            err = 10 * err
        return err
    
    def _get_params(self, model: torch.nn.Module) -> dict:
        """
        Extract current parameter values from a JCAModel.

        Parameters
        ----------
        model : torch.nn.Module
            The JCAModel instance containing learnable parameters.

        Returns
        -------
        params : dict
            Dictionary containing normalized parameter values (fr, phi, tau, vcl, tcl).
        """
        params = {}
        for i, p in enumerate(model.parameters()):
            if i == 0:
                params['fr'] = p.item()
            elif i == 1:
                params['phi'] = p.item()
            elif i == 2:
                params['tau'] = p.item()
            elif i == 3:
                params['vcl'] = p.item()
            elif i == 4:
                params['tcl'] = p.item()
        return params
    
    def _gridsearch(self, base_abs: torch.Tensor, thickness: float) -> tuple:
        """
        Perform a coarse grid search to find good initial parameter estimates for the ML optimizer.

        Parameters
        ----------
        base_abs : torch.Tensor
            Measured absorption coefficients to match.
        thickness : float
            Sample thickness in millimeters.

        Returns
        -------
        tuple
            Normalized initial guesses for (flow_resistivity, porosity, tortuosity, vcl, tcl).
        """
        print("Starting grid search...")
        # Grid search is done on CPU for simplicity (small computation)
        fr = torch.linspace(10000, 1000000, 5)
        phi = torch.linspace(0.05, 0.95, 10)
        tau = torch.linspace(1, 4.5, 5)
        vcl = torch.linspace(10, 450, 10)
        tcl = torch.linspace(10, 450, 10)
        best_err = float('inf')
        best_fr = best_phi = best_tau = best_vcl = best_tcl = None
        
        # Move base_abs to CPU for comparison
        base_abs_cpu = base_abs.cpu() if base_abs.is_cuda else base_abs
        
        for f in fr:
            for p in phi:
                for t in tau:
                    for v in vcl:
                        for tc in tcl:
                            if tc >= v:
                                s = AcousticTMM(incidence='Normal', air_temperature=20, device='cpu')
                                layer = s.Add_JCA_Layer(thickness, f, p, t, v, tc)
                                tm = s.assemble_structure(layer)
                                a = s.absorption(tm)[:, 1].float()
                                err = torch.sum(torch.diff(a - base_abs_cpu) ** 2)
                                
                                if err < best_err:
                                    best_err = err
                                    best_fr = f / 1000000
                                    best_phi = p
                                    best_tau = t / 5
                                    best_vcl = v / 500
                                    best_tcl = tc / 500

        return best_fr, best_phi, best_tau, best_vcl, best_tcl
    
    def ML(self, 
           thickness: float, 
           verbose: bool = True,
           max_iterations: int = 100000,
           learning_rate: float = 1e-3,
           early_stopping_loss: float = 8.0) -> dict:
        """
        Machine learning-based parameter identification using gradient descent optimization.
        
        This method uses the Adam optimizer to find JCA model parameters by minimizing 
        the difference between predicted and measured absorption coefficients. It first 
        performs a coarse grid search to find good initial estimates, then refines them
        using gradient descent.

        Parameters
        ----------
        thickness : float
            The measured thickness of the sample in millimeters.
        verbose : bool, optional
            If True, prints optimization progress every 100 iterations (default is True).
        max_iterations : int, optional
            Maximum number of optimization iterations (default is 100000).
        learning_rate : float, optional
            Initial learning rate for the Adam optimizer (default is 1e-3).
            The learning rate is automatically reduced as the loss decreases.
        early_stopping_loss : float, optional
            Stop optimization when loss falls below this value (default is 8.0).

        Returns
        -------
        result_dict : dict
            Dictionary containing the identified parameters: thickness, flow resistivity,
            porosity, tortuosity, viscous characteristic length, thermal characteristic length,
            and air gap.

        Notes
        -----
        This method requires 'No Gap' mounting condition and uses adaptive learning rate
        scheduling to improve convergence. The learning rate is reduced at loss thresholds
        of 2000, 250, and 10.
        
        Raises
        ------
        ValueError
            If thickness is not positive.
        """
        # Input validation
        if thickness <= 0:
            raise ValueError("Thickness must be positive")
        if max_iterations <= 0:
            raise ValueError("max_iterations must be positive")
        if learning_rate <= 0:
            raise ValueError("learning_rate must be positive")
        
        y = self.meas_abs[0]
        # Move target data to the configured device
        if isinstance(y, torch.Tensor):
            y = y.to(self.device)
        else:
            y = torch.tensor(y, device=self.device)
        
        fr, phi, tau, vcl, tcl1 = self._gridsearch(y, thickness)
        
        # Create model on configured device
        model = JCAModel(fr, phi, tau, vcl, tcl1, self.frequency, device=self.device)
        optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
        
        # Store initial learning rate for adaptive scheduling
        initial_lr = learning_rate

        for t in range(max_iterations):
            y_pred = model.forward(thickness)

            loss = self._criterion(y_pred, y, self._get_params(model))
           
            # Adaptive learning rate scheduling
            if loss < 2000:
                for g in optimizer.param_groups:
                    g['lr'] = initial_lr * 0.5
            if loss < 250:
                for g in optimizer.param_groups:
                    g['lr'] = initial_lr * 0.1
            if loss < 10:
                for g in optimizer.param_groups:
                    g['lr'] = initial_lr * 0.05
            if loss < early_stopping_loss:
                if verbose:
                    print(f"Converged at iteration {t} with loss {loss.item():.4f}")
                break

            if t % 100 == 0 and verbose:
                print(f"Iteration {t}: loss={loss.item():.4f}, {model.string()}")

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

        return model.results()

    def _to_numpy(self, data) -> np.ndarray:
        """
        Convert data to numpy array, handling torch tensors.
        
        Parameters
        ----------
        data : torch.Tensor or np.ndarray or array-like
            Data to convert.
            
        Returns
        -------
        np.ndarray
            Data as numpy array.
        """
        if isinstance(data, torch.Tensor):
            return data.detach().cpu().numpy()
        return np.asarray(data)

    def stats(self, parameters: dict) -> dict:
        """
        Calculates statistics about the parameters identified in the optimization routine via linear regression of the predicted vs measured
        absorption coefficients.

        Parameters
        ----------
        parameters : dict
            Dictionary containing the identified thickness, flow resistivity, porosity, tortuosity, 
            viscous characteristic length, thermal characteristic length, and air gap of the sample.

        Returns
        -------
        stats : dict
            Dictionary containing the slope, intercept, r value, p value, and std error returned from the linear regression.
            If 'Dual' opt_type is specified, the statistics for each mounting condition are averaged.
        """
        if self.opt_type == 'No Gap':
            x = self._to_numpy(self.meas_abs[0])
            y = self._to_numpy(self._predictionJCA(parameters, optimize=False))
            
            slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)
            return {'slope': slope, 'intercept': intercept, 'r_value': r_value, 'p_value': p_value, 'std_err': std_err}
        
        elif self.opt_type == 'Gap':
            x = self._to_numpy(self.meas_abs[1])
            y = self._to_numpy(self._predictionJCA(parameters, optimize=False))
            
            slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)
            return {'slope': slope, 'intercept': intercept, 'r_value': r_value, 'p_value': p_value, 'std_err': std_err}
        
        elif self.opt_type == 'Dual':
            x1 = self._to_numpy(self.meas_abs[0])
            y1 = self._to_numpy(self._predictionJCA(parameters, optimize=False)[0])
            
            slope1, intercept1, r_value1, p_value1, std_err1 = scipy.stats.linregress(x1, y1)
            
            x2 = self._to_numpy(self.meas_abs[1])
            y2 = self._to_numpy(self._predictionJCA(parameters, optimize=False)[1])
            
            slope2, intercept2, r_value2, p_value2, std_err2 = scipy.stats.linregress(x2, y2)
            
            slope = np.abs((slope1 + slope2) / 2)
            intercept = np.abs((intercept1 + intercept2) / 2)
            r_value = np.abs((r_value1 + r_value2) / 2)
            p_value = np.abs((p_value1 + p_value2) / 2)
            std_err = np.abs((std_err1 + std_err2) / 2)
                        
            return {'slope': slope, 'intercept': intercept, 'r_value': r_value, 'p_value': p_value, 'std_err': std_err}

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
            predicted = self._predictionJCA(parameters, optimize=False)
            
            ax.plot(self.frequency,actual,label='Actual')
            ax.plot(self.frequency,predicted,label='Predicted')

            ax.legend(loc="lower right")
            ax.set_ylim(bottom=0)
    
            plt.show()

        
        elif self.opt_type == 'Gap':
            actual = self.meas_abs[1]
            predicted = self._predictionJCA(parameters, optimize=False)
            
            ax.plot(self.frequency,actual,label='Actual')
            ax.plot(self.frequency,predicted,label='Predicted')

            ax.legend(loc="lower right")
            ax.set_ylim(bottom=0)
    
            plt.show()
        
        elif self.opt_type == 'Dual':
            no_gap_pred = self._predictionJCA(parameters, optimize=False)[0]
            gap_pred = self._predictionJCA(parameters, optimize=False)[1]
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
            predicted = self._predictionJCA(parameters, optimize=False)
            all_data = {'frequency':self.frequency,'measured':actual,'predicted':predicted}

        
        elif self.opt_type == 'Gap':
            actual = self.meas_abs[1]
            predicted = self._predictionJCA(parameters, optimize=False)
            all_data = {'frequency':self.frequency,'measured':actual,'predicted':predicted}
        
        elif self.opt_type == 'Dual':
            no_gap_pred = self._predictionJCA(parameters, optimize=False)[0]
            gap_pred = self._predictionJCA(parameters, optimize=False)[1]
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

class JCAModel(AcousticTMM):
    """
    PyTorch model for JCA parameter optimization using gradient descent.
    
    This model wraps AcousticTMM to enable gradient-based optimization of
    JCA material parameters.
    
    Parameters
    ----------
    best_fr : torch.Tensor
        Initial normalized flow resistivity guess.
    best_phi : torch.Tensor
        Initial porosity guess.
    best_tau : torch.Tensor
        Initial normalized tortuosity guess.
    best_vcl : torch.Tensor
        Initial normalized viscous characteristic length guess.
    best_tcl : torch.Tensor
        Initial normalized thermal characteristic length guess.
    freq : torch.Tensor
        Frequency array for calculations.
    device : str, optional
        Device for computations ('cpu' or 'cuda'). Default is 'cpu'.
    """
    
    def __init__(self,
                 best_fr,
                 best_phi,
                 best_tau,
                 best_vcl,
                 best_tcl,
                 freq,
                 device: str = 'cpu'):

        super().__init__(device=device)
        self._device = device
        
        # Ensure initial values are tensors on the correct device
        self.fr = torch.nn.Parameter(torch.tensor(best_fr, device=device) if not isinstance(best_fr, torch.Tensor) else best_fr.to(device))
        self.phi = torch.nn.Parameter(torch.tensor(best_phi, device=device) if not isinstance(best_phi, torch.Tensor) else best_phi.to(device))
        self.tau = torch.nn.Parameter(torch.tensor(best_tau, device=device) if not isinstance(best_tau, torch.Tensor) else best_tau.to(device))
        self.vcl = torch.nn.Parameter(torch.tensor(best_vcl, device=device) if not isinstance(best_vcl, torch.Tensor) else best_vcl.to(device))
        self.tcl = torch.nn.Parameter(torch.tensor(best_tcl, device=device) if not isinstance(best_tcl, torch.Tensor) else best_tcl.to(device))
        
        self.structure = AcousticTMM(incidence='Normal', air_temperature=20, device=device)
        self.structure.frequency = freq
        self.thickness = None
    
    def forward(self, thickness: float) -> torch.Tensor:
        """
        Compute absorption coefficients for current parameters.
        
        Parameters
        ----------
        thickness : float
            Sample thickness in millimeters.
            
        Returns
        -------
        torch.Tensor
            Absorption coefficients at each frequency.
        """
        self.thickness = thickness
        layer = self.structure.Add_JCA_Layer(
            thickness, 
            self.fr * 1000000, 
            self.phi, 
            self.tau * 5, 
            self.vcl * 500, 
            self.tcl * 500
        )
        tm = self.structure.assemble_structure(layer)
        a = self.structure.absorption(tm)[:, 1].float()
        return a
    
    def string(self) -> str:
        """Return string representation of current parameters."""
        return (f'fr = {self.fr.item()*1000000:.1f} '
                f'phi = {self.phi.item():.4f} '
                f'tau = {self.tau.item()*5:.3f} '
                f'vcl = {self.vcl.item()*500:.1f} '
                f'tcl = {self.tcl.item()*500:.1f}')
    
    def results(self) -> dict:
        """Return identified parameters as a dictionary."""
        return {
            'thickness': self.thickness,
            'flow resistivity': self.fr.item() * 1000000,
            'porosity': self.phi.item(),
            'tortuosity': self.tau.item() * 5,
            'viscous characteristic length': self.vcl.item() * 500,
            'thermal characteristic length': self.tcl.item() * 500,
            'air gap': 0
        }