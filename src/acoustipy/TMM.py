import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from acoustipy.Database import AcoustiBase


class AcousticTMM():
    """
    Create an AcousticTMM object
    
    Description:
    ------------
    AcousticTMM is an implementation of the acoustic transfer matrix method and a number of porous material models built on top of numpy (https://numpy.org/).  
    It can be used to calculate interesting acoustic characteristics -- like the frequency dependent reflection, absorption, and transmission coefficients of a variety of materials. 
    Using AcousticTMM, a number of layers can be defined and combined to create multilayer structures, which can then be acoustically simulated via the 
    transfer matrix method.
    
    Attributes:
    -----------
    fmin (int):
       Minimum frequency of the range of interest

    fmax (int)::
        Maximum frequency of the range of interest

    fs (int):
        Frequency step size; reduce computation time by taking larger strides between frequencies within the range.

    incidence (str):
        'Normal' or 'Diffuse'; Compute properties of interest in impedence tube like conditions or in reverberant conditions.

    angles (list):
        If Diffuse incidence is specified, gives tighter control over the angles used to calculate properties of interest.

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
                 fmin:int = 4,
                 fmax:int = 3200,
                 fs:int = 4,
                 incidence:str = 'Normal',
                 angles:list[int] = [0,79,1],
                 air_temperature:float = None,
                 sound_speed:float = 343.152,
                 air_density:float = 1.2058,
                 Cp:float = 1.004425,
                 Cv:float = 0.717425,
                 viscosity:float = 1.825e-05,
                 Pr:float = .7157,
                 P0:float = 101325
                 ):
        
        #Define limits on the frequency range and angles
        if fs <= 0:
            raise ValueError("Frequency Step (fs) must be greater than 0!")
        
        if fmin <= 0:
            raise ValueError("Minimum frequency (fmin) must be greater than 0!")
        
        if fmin >= fmax:
            raise ValueError("Minimum frequency (fmin) must be less than maximum frequency (fmax)!")
            
        if angles[0] >= angles[1]:
            raise ValueError("Minimum angle must be less than maximum angle!")
        
        if abs(round((angles[1]-angles[0])%angles[2]/angles[2])-((angles[1]-angles[0])%angles[2]/angles[2])) >  1e-6:
            raise ValueError("Angle step size must be a multiple of the angle range!")
        
        
        self.temp = air_temperature
        self.speed = sound_speed
        self.density = air_density
        self.fmin = fmin
        self.fmax = fmax
        self.fs = fs
        self.incidence = incidence
        self.angles = angles
        self.THIRD_OCTAVE_PREFERRED = np.asarray([16,20,25,31.5,40,50,63,
                                                  80,100,125,160,200,250,
                                                  315,400,500,630,800,1000,
                                                  1250,1600,2000,2500,3150,
                                                  4000,5000,6300,8000,10000,
                                                  12500,16000,20000])
        self.OCTAVE_PREFERRED  = self.THIRD_OCTAVE_PREFERRED[0::3]
        self.Cp = Cp
        self.Cv = Cv
        self.viscosity = viscosity
        self.Pr = Pr
        self.P0 = P0
        self._custom_freq =  np.arange(self.fmin,self.fmax+self.fs,self.fs)
        self.layers = []
        

    @property
    def frequency(self):
        #frequency range of interest
        return(self._custom_freq)
    
    @frequency.setter
    def frequency(self, value):
        self._custom_freq = value

    @property
    def ang_freq(self):
        
        #angular frequency
        w = 2.0*np.pi*self.frequency
        return(w)
        
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
    def k0(self):
        #Wavenumber
        k0 = self.ang_freq / self.soundspeed_temp
        return (k0)
    
    @property
    def Z0(self):
        #Characteristic Impedance of Air
        z0 = self.density_temp*self.soundspeed_temp
        return(z0)

    
    def _create_layer_TM(self,
                         Zp: np.ndarray,
                         kp: np.ndarray,
                         thickness: float) -> np.ndarray:
        """
        Creates the transfer matrix for an individual layer in either normal or diffuse sound fields

        Parameters
        ----------
        Zp (ndarray):
            Characteristic impedance of the layer
            
        kp (ndarray):
            Characteristic wavenumber of the layer

        thickness (float):
            Layer thickness [m]

        Returns
        -------
        TM (ndarray):
            Normal incidence --> 2 x 2 x len(frequency) numpy array representing the transfer matrix
            Diffuse field --> 2 x 2 x len(frequency) x len(angles) array representing the transfer matrix
        
        """
        if self.incidence == "Normal":
            TM = np.zeros((2,2,len(self.frequency)),dtype = 'complex_')
            TM[0][0] = np.cos(kp*thickness)
            TM[0][1] = 1j*Zp*np.sin(kp*thickness)
            TM[1][0] = (1j/Zp)*np.sin(kp*thickness)
            TM[1][1] = np.cos(kp*thickness)
            
            return (TM)

        elif self.incidence == "Diffuse":
            angles = np.arange(self.angles[0],self.angles[1],self.angles[2])
            vs = np.sin(np.radians(angles))
            TM = np.zeros((2,2,len(self.frequency),len(angles)),dtype = 'complex_')
            kpx = np.zeros((len(self.frequency),len(angles)),dtype = 'complex_')
            
            Kp = np.einsum('ij,ij -> ij',np.tile(kp[:,None], (1,len(angles))),np.tile(kp[:,None], (1,len(angles))))
            K0 = np.einsum('ij,ij -> ij',np.tile(self.k0[:,None],(1,len(angles))),np.tile(self.k0[:,None],(1,len(angles))))
            VS = np.einsum('ij,ij -> ij',vs[None,:],vs[None,:])
            kpx[:,:] = np.sqrt((Kp)-(np.einsum('ij,ij -> ij',K0,VS)))

            sin = 1j*np.sin(kpx*thickness)
            cos = np.cos(kpx*thickness)
            offset = np.einsum('ij,ij,ij -> ij',Zp[:, None],kp[:, None],1/kpx)

            
            TM[0,0,:,:] = cos
            TM[0,1,:,:] = np.einsum('ij,ij -> ij',offset,sin)
            TM[1,0,:,:] = np.einsum('ij,ij -> ij',sin,1/offset)
            TM[1,1,:,:] = cos
            
            return(TM)
        
    def _create_Maa_MPP_TM(self,
                           Zp: np.ndarray) -> np.ndarray:
        """
        Creates the transfer matrix for a Maa microperforated panel layer in either normal or diffuse sound fields

        Parameters
        ----------
        Zp (ndarray):
            Characteristic impedance of the layer

        Returns
        -------
        TM (ndarray):
            Normal incidence --> 2 x 2 x len(frequency) numpy array representing the transfer matrix
            Diffuse field --> 2 x 2 x len(frequency) x len(angles) array representing the transfer matrix
        
        """
        if self.incidence == "Normal":
            TM = np.zeros((2,2,len(self.frequency)),dtype = 'complex_')
            TM[0][0] = 1
            TM[0][1] = Zp
            TM[1][0] = 0
            TM[1][1] = 1
            
            return (TM)
        
        elif self.incidence == "Diffuse":
            angles = np.arange(self.angles[0],self.angles[1],self.angles[2])
            
            TM = np.zeros((2,2,len(self.frequency),len(angles)),dtype = 'complex_')

            count = 0
            for theta in angles:
               
                TM[0,0,:,count] = 1
                TM[0,1,:,count] = Zp*np.cos(np.radians(theta))
                TM[1,0,:,count] = 0
                TM[1,1,:,count] = 1
                count += 1

            return(TM)
    
    def _calc_dynamics(self,
                       flow_resistivity: float,
                       porosity: float,
                       tortuosity: float,
                       viscous_characteristic_length: float,
                       thermal_characteristic_length: float,
                       thermal_permeability: float,
                       thermal_tortuosity: float,
                       viscous_tortuosity: float,
                       model: str) -> tuple[np.ndarray, np.ndarray]:
        """
        Calculates the dynamic mass density and bulk modulus for porous, equivalent fluid models.

        Parameters
        ----------
        flow_resisitivty (float):
            Static air flow resistivity of the material [Pa*s/m2]
            
        porosity (float):
            open porosity of the material

        tortuosity (float):
            high frequency limit of the tortuosity of the material

        viscous_characteristic_length (float):
            viscous characteristic length of the material [m]

        thermal_characteristic_length (float):
            thermal characteristic length of the material [m]

        thermal_permeability (float):
            static thermal permeability of the material [m2]

        thermal_tortuosity (float):
            static thermal tortuosity of the material

        viscous_tortuosity (float):
            static viscous tortuosity of the material

        model (str):
            defines the equivalent fluid model to be used
        

        Returns
        -------
        peff, keff (tuple(ndarray, ndarray)):
            arrays of shape [len(frequency),1] representing the dynamic mass density
            and bulk modulus of the material

        """
        fr = flow_resistivity
        tau = tortuosity
        phi = porosity
        vcl = viscous_characteristic_length
        tcl = thermal_characteristic_length
        kprime = thermal_permeability
        tauprime = thermal_tortuosity
        tau0 = viscous_tortuosity
        
        w = self.ang_freq
        
        if model == 'DB':
            #Delaney-Bazley Model
            DB_Param = (self.frequency*self.density_temp)/fr
        
            a = DB_Param**-.754
            b = DB_Param**-.732
            c = DB_Param**-.700
            d = DB_Param**-.595
            
            Zr = self.Z0*(1+.0571*a)
            Zi = self.Z0*.0870*b
            kr = self.k0*(1+.0978*c)
            ki = .189*self.k0*d
            
            
            Zp = Zr-(1j*Zi)
            kp = kr-(1j*ki)

            peff = (Zp*kp)/w
            keff = (Zp*w)/kp

        elif model == 'DBM':
            #Delaney-Bazley-Miki Model
            DBM_Param = self.frequency/fr
        
            a = DBM_Param**-.632
            b = DBM_Param**-.618
            
            Zr = self.Z0*(1+.0699*a)
            Zi = self.Z0*.107*a
            kr = self.k0*(1+.109*b)
            ki = .160*self.k0*b

            Zp = Zr-(1j*Zi)
            kp = kr-(1j*ki)

            peff = (Zp*kp)/w
            keff = (Zp*w)/kp
        
        elif model == 'JCA':
            #Johnson-Champoux-Allard Model
            wp = (w*self.density_temp*tau)/(fr*phi)
            wpp = ((tcl**2)*self.Pr_temp*self.density_temp*w)/(8*self.viscosity_temp)
            p=1
            m = (8*self.viscosity_temp*tau)/(fr*phi*(vcl**2))
            pprime = 1
            mprime = 1
            
            Fw = 1-p+(p*np.sqrt(1+((1j*m*wp)/(2*(p**2)))))
            aw = tau*(1+(Fw/(1j*wp)))
            
            Fwp = 1-pprime+(pprime*np.sqrt(1+((1j*mprime*wpp)/(2*(pprime**2)))))
            bw = self.gamma_temp-((self.gamma_temp-1)*((1+(Fwp/(1j*wpp)))**-1)) 
            
            peff = self.density_temp*aw/phi
            keff = ((self.gamma_temp*self.P0)/(phi*bw))
        
        elif model == 'JCAL':
            #Johnson-Champoux-Allard-Lafarge Model
            wp = (w*self.density_temp*tau)/(fr*phi)
            wpp = (w*self.density_temp*self.Pr_temp*kprime)/(self.viscosity_temp*phi)
            p = 1
            m = (8*self.viscosity_temp*tau)/(fr*phi*(vcl**2))
            pprime = 1
            mprime = (8*kprime)/(phi*(tcl**2))

            Fw = 1-p+(p*np.sqrt(1+((1j*m*wp)/(2*(p**2)))))
            aw = tau*(1+(Fw/(1j*wp)))
            
            Fwp = 1-pprime+(pprime*np.sqrt(1+((1j*mprime*wpp)/(2*(pprime**2)))))
            bw = self.gamma_temp-((self.gamma_temp-1)*((1+(Fwp/(1j*wpp)))**-1)) 
            
            peff = self.density_temp*aw/phi
            keff = ((self.gamma_temp*self.P0)/(phi*bw))
            
        elif model == 'JCAPL':
            #Johnson-Champoux-Allard-Pride-Lafarge Model
            wp = (w*self.density_temp*tau)/(fr*phi)
            wpp = (w*self.density_temp*self.Pr_temp*kprime)/(self.viscosity_temp*phi)
            m = (8*self.viscosity_temp*tau)/(fr*phi*(vcl**2))
            p = m/(4*((tau0/tau)-1))
            mprime = (8*kprime)/(phi*(tcl**2))
            pprime = mprime/((4*(tauprime-1)))
            
            Fw = 1-p+(p*np.sqrt(1+((1j*m*wp)/(2*(p**2)))))
            aw = tau*(1+(Fw/(1j*wp)))
            
            Fwp = 1-pprime+(pprime*np.sqrt(1+((1j*mprime*wpp)/(2*(pprime**2)))))
            bw = self.gamma_temp-((self.gamma_temp-1)*((1+(Fwp/(1j*wpp)))**-1)) 
            
            peff = self.density_temp*aw/phi
            keff = ((self.gamma_temp*self.P0)/(phi*bw))
        
        return (peff,keff)

    def Add_Air_Layer(self,
                      thickness: float=400,
                      save_layer: bool = False,
                      layer_name: str = None) -> list[np.ndarray, float, str]:
        """
        Define an air gap layer

        Parameters
        ----------
        thickness (float):
            air gap thickness [mm]

        save_layer (bool):
            Specify whether to save the input parameters to a database for later use.

        layer_name (str):
            If save_layer is set to True, specify the name of the layer.  Must be a unique identifier.

        Returns
        -------
        TM, thickness, layer_name (list[ndarray, float, str]):
            The transfer matrix, thickness, and name of the layer.
        
        """
        if save_layer == True:
            params = [layer_name,'AIR','null',thickness,'null','null','null','null','null','null','null','null','null','null','null','null','null']
            self._layer_to_db(params)

        Zp = np.full(len(self.frequency),self.Z0)
        thickness = thickness / 1000
        TM = self._create_layer_TM(Zp,self.k0,thickness)
        
        return([TM,0,layer_name])
    
    def Add_DB_Layer(self,
                     thickness: float,
                     flow_resistivity: float,
                     save_layer: bool = False,
                     layer_name: str = None) -> list[np.ndarray, float, str]:
        """
        Define a layer using the Delaney-Bazley Model

        Parameters
        ----------
        thickness (float):
            layer thickness [mm]
        
        flow_resistivity (float):
            static flow resistivity of the layer [Pa*s/m2]
        
        save_layer (bool):
            Specify whether to save the input parameters to a database for later use.

        layer_name (str):
            If save_layer is set to True, specify the name of the layer.  Must be a unique identifier.

        Returns
        -------
        TM, thickness, layer_name (list(ndarray, float, str):
            The transfer matrix, thickness, and name of the layer.
        
        """      
        if save_layer == True:
            params = [layer_name,'DB','null',thickness,flow_resistivity,'null','null','null','null','null','null','null','null','null','null','null','null']
            self._layer_to_db(params)
        
        thickness = thickness/1000
        
        dyns = self._calc_dynamics(flow_resistivity, 0, 0, 0, 0, 0, 0, 0, model='DB')
      
        Zp = np.sqrt(dyns[0]*dyns[1])
        kp = self.ang_freq*np.sqrt(dyns[0]/dyns[1])

        TM = self._create_layer_TM(Zp,kp,thickness)

        return([TM,thickness,layer_name])
    
    def Add_DBM_Layer(self,
                      thickness: float,
                      flow_resistivity: float,
                      save_layer: bool = False,
                      layer_name: str = None) -> list[np.ndarray, float, str]:
        """
        Define a layer using the Delaney-Bazley-Miki Model

        Parameters
        ----------
        thickness (float):
            layer thickness [mm]
        
        flow_resistivity (float):
            static flow resistivity of the layer [Pa*s/m2]

        save_layer (bool):
            Specify whether to save the input parameters to a database for later use.

        layer_name (str):
            If save_layer is set to True, specify the name of the layer.  Must be a unique identifier.

        Returns
        -------
        TM, thickness, layer_name (list[ndarray, float, str]):
            The transfer matrix, thickness, and name of the layer.
        
        """
        if save_layer == True:
            params = [layer_name,'DBM','null',thickness,flow_resistivity,'null','null','null','null','null','null','null','null','null','null','null','null']
            self._layer_to_db(params)
            
        thickness = thickness/1000

        dyns = self._calc_dynamics(flow_resistivity, 0, 0, 0, 0, 0, 0, 0, model='DB')
      
        Zp = np.sqrt(dyns[0]*dyns[1])
        kp = self.ang_freq*np.sqrt(dyns[0]/dyns[1])
        
        TM = self._create_layer_TM(Zp,kp,thickness)
        
        return([TM,thickness,layer_name])
       
    
    def Add_JCA_Layer(self,
                      thickness: float,
                      flow_resistivity: float,
                      porosity: float,
                      tortuosity: float,
                      viscous_characteristic_length: float,
                      thermal_characteristic_length: float,
                      save_layer: bool = False,
                      layer_name: str = None) -> list[np.ndarray, float, str]:
        """
        Define a layer using the Johnson-Champoux-Allard Model

        Parameters
        ----------
        thickness (float):
            layer thickness [mm]
        
        flow_resisitivty (float):
            Static air flow resistivity of the material [Pa*s/m2]
            
        porosity (float):
            open porosity of the material

        tortuosity (float):
            high frequency limit of the tortuosity of the material

        viscous_characteristic_length (float):
            viscous characteristic length of the material [µm]

        thermal_characteristic_length (float):
            thermal characteristic length of the material [µm]

        save_layer (bool):
            Specify whether to save the input parameters to a database for later use.

        layer_name (str):
            If save_layer is set to True, specify the name of the layer.  Must be a unique identifier.
        
        Returns
        -------
        TM, thickness, layer_name (list[ndarray, float, str]):
            The transfer matrix, thickness, and name of the layer.
        
        """
        if save_layer == True:
            params = [layer_name,'JCA','null',thickness,flow_resistivity,porosity,tortuosity,viscous_characteristic_length,thermal_characteristic_length,'null','null','null','null','null','null','null','null']
            self._layer_to_db(params)
        
        thickness = thickness/1000
        
        fr = flow_resistivity
        tau = tortuosity
        phi = porosity
        vcl = viscous_characteristic_length*(1e-6)
        tcl = thermal_characteristic_length*(1e-6)
        
        dyns = self._calc_dynamics(fr, phi, tau, vcl, tcl, 0, 0, 0, model='JCA')
      
        Zp = np.sqrt(dyns[0]*dyns[1])
        kp = self.ang_freq*np.sqrt(dyns[0]/dyns[1])
        
        TM = self._create_layer_TM(Zp,kp,thickness)
         
        return([TM,thickness,layer_name])
    
    def Add_JCAL_Layer(self,
                       thickness: float,
                       flow_resistivity: float,
                       porosity: float,
                       tortuosity: float,
                       viscous_characteristic_length: float,
                       thermal_characteristic_length: float,
                       thermal_permeability: float,
                       save_layer: bool = False,
                       layer_name: str = None) -> list[np.ndarray, float, str]:
        """
        Define a layer using the Johnson-Champoux-Allard-Lafarge Model

        Parameters
        ----------
        thickness (float):
            layer thickness [mm]
        
        flow_resisitivty (float):
            Static air flow resistivity of the material [Pa*s/m2]
            
        porosity (float):
            open porosity of the material

        tortuosity (float):
            high frequency limit of the tortuosity of the material

        viscous_characteristic_length (float):
            viscous characteristic length of the material [µm]

        thermal_characteristic_length (float):
            thermal characteristic length of the material [µm]

        thermal_permeability (float):
            static thermal permeability of the material [m2]

        save_layer (bool):
            Specify whether to save the input parameters to a database for later use.

        layer_name (str):
            If save_layer is set to True, specify the name of the layer.  Must be a unique identifier.
        
        Returns
        -------
        TM, thickness, layer_name (list[ndarray, float, str]):
            The transfer matrix, thickness, and name of the layer.
        
        """
        if save_layer == True:
            params = [layer_name,'JCAL','null',thickness,flow_resistivity,porosity,tortuosity,viscous_characteristic_length,thermal_characteristic_length,thermal_permeability,'null','null','null','null','null','null','null']
            self._layer_to_db(params)
            
        thickness = thickness/1000
        
        fr = flow_resistivity
        tau = tortuosity
        phi = porosity
        vcl = viscous_characteristic_length*(1e-6)
        tcl = thermal_characteristic_length*(1e-6)
        kprime = thermal_permeability*(1e-10)

        dyns = self._calc_dynamics(fr, phi, tau, vcl, tcl, kprime, 0, 0, model='JCAL')
        
        Zp = np.sqrt(dyns[0]*dyns[1])
        kp = self.ang_freq*np.sqrt(dyns[0]/dyns[1])
        
        TM = self._create_layer_TM(Zp,kp,thickness)
        
        return([TM,thickness,layer_name])
    
    def Add_JCAPL_Layer(self,
                        thickness: float,
                        flow_resistivity: float,
                        porosity: float,
                        tortuosity: float,
                        viscous_characteristic_length: float,
                        thermal_characteristic_length: float,
                        thermal_permeability: float,
                        thermal_tortuosity: float,
                        viscous_tortuosity: float,
                        save_layer: bool = False,
                        layer_name: str = None) -> list[np.ndarray, float, str]:
        """
        Define a layer using the Johnson-Champoux-Allard-Pride-Lafarge Model

        Parameters
        ----------
        thickness (float):
            layer thickness [mm]
        
        flow_resisitivty (float):
            Static air flow resistivity of the material [Pa*s/m2]
            
        porosity (float):
            open porosity of the material

        tortuosity (float):
            high frequency limit of the tortuosity of the material

        viscous_characteristic_length (float):
            viscous characteristic length of the material [µm]

        thermal_characteristic_length (float):
            thermal characteristic length of the material [µm]

        thermal_permeability (float):
            static thermal permeability of the material [m2]

        thermal_tortuosity (float):
            static thermal tortuosity of the material

        viscous_tortuosity (float):
            static viscous tortuosity of the material

        save_layer (bool):
            Specify whether to save the input parameters to a database for later use.

        layer_name (str):
            If save_layer is set to True, specify the name of the layer.  Must be a unique identifier.
        
        Returns
        -------
        TM, thickness, layer_name (list[ndarray, float, str]):
            The transfer matrix, thickness, and name of the layer.
        
        """
        if save_layer == True:
            params = [layer_name,'JCAPL','null',thickness,flow_resistivity,porosity,tortuosity,viscous_characteristic_length,thermal_characteristic_length,thermal_permeability,thermal_tortuosity,viscous_tortuosity,'null','null','null','null','null']
            self._layer_to_db(params)
        
        thickness = thickness/1000

        fr = self.viscosity_temp/flow_resistivity
        tau = tortuosity
        phi = porosity
        vcl = viscous_characteristic_length*(1e-6)
        tcl = thermal_characteristic_length*(1e-6)
        kprime = thermal_permeability*(1e-10)
        tauprime = thermal_tortuosity
        tau0 = viscous_tortuosity
        
        dyns = self._calc_dynamics(fr, phi, tau, vcl, tcl, kprime, tauprime, tau0, model='JCAPL')
        
        Zp = np.sqrt(dyns[0]*dyns[1])
        kp = self.ang_freq*np.sqrt(dyns[0]/dyns[1])

        TM = self._create_layer_TM(Zp,kp,thickness)
        
        return([TM,thickness,layer_name])
    
    def Add_Horoshenkov_Layer(self,
                              thickness: float,
                              porosity: float,
                              median_pore_size: float,
                              pore_size_distribution: float,
                              save_layer: bool = False,
                              layer_name: str = None) -> list[np.ndarray, float, str]:
        """
        Define a layer using the Horoshenkov et al Model

        https://pubs.aip.org/asa/jasa/article/145/4/2512/845598/A-three-parameter-analytical-model-for-the
        
        Parameters
        ----------
        thickness (float):
            layer thickness [mm]
            
        porosity (float):
            open porosity of the material

        median_pore_size (float):
            median pore size of the material [µm]
        
        pore_size_distribution (float):
            standard deviation in the pore size distribution

        save_layer (bool):
            Specify whether to save the input parameters to a database for later use.

        layer_name (str):
            If save_layer is set to True, specify the name of the layer.  Must be a unique identifier.
        
        Returns
        -------
        TM, thickness, layer_name (list[ndarray, float, str]):
            The transfer matrix, thickness, and name of the layer.
        
        """
        if save_layer == True:
            params = [layer_name,'horoshenkov','null',thickness,'null',porosity,'null','null','null','null','null','null','null','null','null',median_pore_size,pore_size_distribution]
            self._layer_to_db(params)
            
        thickness = thickness/1000
        mps = median_pore_size*(1e-6)
        psd = pore_size_distribution
        phi = porosity

        X = (psd*np.log(2))**2
        tau = np.exp(4*X)
        Z = phi*(mps**2)/(8*tau)
        
        
        fr = (self.viscosity_temp/Z)*np.exp(6*X)
        vcl = mps*np.exp((-5/2)*X)
        tcl = mps*np.exp((3/2)*X)
        kprime = Z/np.exp(-6*X)
        
        dyns = self._calc_dynamics(fr, phi, tau, vcl, tcl, kprime, 0, 0, model='JCAL')
        
        Zp = np.sqrt(dyns[0]*dyns[1])
        kp = self.ang_freq*np.sqrt(dyns[0]/dyns[1])
        
        TM = self._create_layer_TM(Zp,kp,thickness)
        
        return([TM,thickness,layer_name])
    
    def Add_Biot_Limp_Layer(self,
                            EF_model: str,
                            thickness: float,
                            flow_resistivity: float,
                            mass_density: float,
                            porosity: float,
                            tortuosity: float=0,
                            viscous_characteristic_length: float=0,
                            thermal_characteristic_length: float=0,
                            thermal_permeability: float=0,
                            thermal_tortuosity: float=0,
                            viscous_tortuosity: float=0,
                            save_layer: bool = False,
                            layer_name: str = None) -> list[np.ndarray, float, str]:
                 
        """
        Define a limp Biot layer, using any of the equivalent fluid models
        
        https://doi.org/10.1121/1.4826175 F.-X. Bécot, L. Jaouen, An alternative Biot's formulation for dissipative porous media with skeleton deformation, J. Acoust. Soc. Am. 134 (6), pp. 4801-4807, 2013

        Parameters
        ----------
        EF_model (str):
            Equivalent fluid model to be used:
            DB --> Delaney-Bazley
            DBM --> Delaney-Bazley-Miki
            JCA --> Johnson-Champoux-Allard
            JCAL --> Johnson-Champoux-Allar-Lafarge
            JCAPL --> Johnson-Champoux-Allard-Pride-Lafarge

        thickness (float):
            layer thickness [mm]
        
        flow_resisitivty (float):
            Static air flow resistivity of the material [Pa*s/m2]
            
        porosity (float):
            open porosity of the material

        mass_density (float):
            bulk density of the material [kg/m3]

        tortuosity (float):, optional
            high frequency limit of the tortuosity of the material. Needed for JCA, JCAL, and JCAPL models.

        viscous_characteristic_length (float):, optional
            viscous characteristic length of the material [µm]. Needed for JCA, JCAL, and JCAPL models.

        thermal_characteristic_length (float):, optional
            thermal characteristic length of the material [µm]. Needed for JCA, JCAL and JCAPL models.

        thermal_permeability (float):, optional
            static thermal permeability of the material [m2]. Needed for JCAL and JCAPL models.

        thermal_tortuosity (float):, optional
            static thermal tortuosity of the material. Needed for JCAPL models.

        viscous_tortuosity (float):, optional
            static viscous tortuosity of the material. Needed for JCAPL models.

        save_layer (bool):
            Specify whether to save the input parameters to a database for later use.

        layer_name (str):
            If save_layer is set to True, specify the name of the layer.  Must be a unique identifier.
        
        Returns
        -------
        TM, thickness, layer_name (list[ndarray, float, str]):
            The transfer matrix, thickness, and name of the layer.
        
        """
        if save_layer == True:
            params = [layer_name,'biot_limp',EF_model,thickness,flow_resistivity,porosity,tortuosity,viscous_characteristic_length,thermal_characteristic_length,thermal_permeability,thermal_tortuosity,viscous_tortuosity,mass_density,'null','null','null','null']
            self._layer_to_db(params)
        
        thickness = thickness/1000

        fr = flow_resistivity
        tau = tortuosity
        phi = porosity
        vcl = viscous_characteristic_length*(1e-6)
        tcl = thermal_characteristic_length*(1e-6)
        kprime = thermal_permeability*(1e-10)
        tauprime = thermal_tortuosity
        tau0 = viscous_tortuosity
        
        if EF_model == 'JCA':
            peff,keff = self._calc_dynamics(fr, phi, tau, vcl, tcl, kprime, tauprime, tau0, model=EF_model)
        
        elif EF_model == 'JCAL':
            peff,keff = self._calc_dynamics(fr, phi, tau, vcl, tcl, kprime, tauprime, tau0, model=EF_model)
        
        elif EF_model == 'JCAPL':
            fr = self.viscosity_temp/fr
            peff,keff = self._calc_dynamics(fr, phi, tau, vcl, tcl, kprime, tauprime, tau0, model=EF_model)
        
        elif EF_model == 'DB':
            peff,keff = self._calc_dynamics(fr, phi, tau, vcl, tcl, kprime, tauprime, tau0, model=EF_model)

        elif EF_model == 'DBM':
            peff,keff = self._calc_dynamics(fr, phi, tau, vcl, tcl, kprime, tauprime, tau0, model=EF_model)

        rho_tilde = mass_density+(phi*self.density_temp)-((self.density_temp**2)/peff)       ##eq. 15
        
        gamma_tilde = (self.density_temp/peff)-1     #eq. 16
        
        rho_eq_limp = (1/(phi*peff)) + ((gamma_tilde**2)/(phi*rho_tilde))    #eq. 24
        
        rho_eq_limp = 1/rho_eq_limp
        
        Zp = np.sqrt(rho_eq_limp*keff)
        kp = self.ang_freq*np.sqrt(rho_eq_limp/keff)

        TM = self._create_layer_TM(Zp,kp,thickness)
        
        return([TM,thickness,layer_name])
    
    def Add_Biot_Rigid_Layer(self,
                             EF_model: str,
                             thickness: float,
                             flow_resistivity: float,
                             mass_density: float,
                             porosity: float,
                             tortuosity: float=0,
                             viscous_characteristic_length: float=0,
                             thermal_characteristic_length: float=0,
                             thermal_permeability: float=0,
                             thermal_tortuosity: float=0,
                             viscous_tortuosity: float=0,
                             save_layer: bool = False,
                             layer_name: str = None) -> list[np.ndarray, float, str]:
                
        """
        Define a rigid Biot layer, using any of the equivalent fluid models

        https://doi.org/10.1121/1.4826175 F.-X. Bécot, L. Jaouen, An alternative Biot's formulation for dissipative porous media with skeleton deformation, J. Acoust. Soc. Am. 134 (6), pp. 4801-4807, 2013


        Parameters
        ----------
        EF_model (str):
            Equivalent fluid model to be used:
            DB --> Delaney-Bazley
            DBM --> Delaney-Bazley-Miki
            JCA --> Johnson-Champoux-Allard
            JCAL --> Johnson-Champoux-Allar-Lafarge
            JCAPL --> Johnson-Champoux-Allard-Pride-Lafarge

        thickness (float):
            layer thickness [mm]
        
        flow_resisitivty (float):
            Static air flow resistivity of the material [Pa*s/m2]
            
        porosity (float):
            open porosity of the material

        mass_density (float):
            bulk density of the material [kg/m3]

        tortuosity (float):, optional
            high frequency limit of the tortuosity of the material. Needed for JCA, JCAL, and JCAPL models.

        viscous_characteristic_length (float):, optional
            viscous characteristic length of the material [µm]. Needed for JCA, JCAL, and JCAPL models.

        thermal_characteristic_length (float):, optional
            thermal characteristic length of the material [µm]. Needed for JCA, JCAL and JCAPL models.

        thermal_permeability (float):, optional
            static thermal permeability of the material [m2]. Needed for JCAL and JCAPL models.

        thermal_tortuosity (float):, optional
            static thermal tortuosity of the material. Needed for JCAPL models.

        viscous_tortuosity (float):, optional
            static viscous tortuosity of the material. Needed for JCAPL models.

        save_layer (bool):
            Specify whether to save the input parameters to a database for later use.

        layer_name (str):
            If save_layer is set to True, specify the name of the layer.  Must be a unique identifier.
        
        Returns
        -------
        TM, thickness, layer_name (list[ndarray, float, str]):
            The transfer matrix, thickness, and name of the layer.
        
        """
        if save_layer == True:
            params = [layer_name,'biot_rigid',EF_model,thickness,flow_resistivity,porosity,tortuosity,viscous_characteristic_length,thermal_characteristic_length,thermal_permeability,thermal_tortuosity,viscous_tortuosity,mass_density,'null','null','null','null']
            self._layer_to_db(params)
        
        thickness = thickness/1000

        fr = flow_resistivity
        tau = tortuosity
        phi = porosity
        vcl = viscous_characteristic_length*(1e-6)
        tcl = thermal_characteristic_length*(1e-6)
        kprime = thermal_permeability*(1e-10)
        tauprime = thermal_tortuosity
        tau0 = viscous_tortuosity

        
        if EF_model == 'JCA':
            peff,keff = self._calc_dynamics(fr, phi, tau, vcl, tcl, kprime, tauprime, tau0, model=EF_model)
        
        elif EF_model == 'JCAL':
            peff,keff = self._calc_dynamics(fr, phi, tau, vcl, tcl, kprime, tauprime, tau0, model=EF_model)
        
        elif EF_model == 'JCAPL':
            fr = self.viscosity_temp/fr
            peff,keff = self._calc_dynamics(fr, phi, tau, vcl, tcl, kprime, tauprime, tau0, model=EF_model)
        
        elif EF_model == 'DB':
            peff,keff = self._calc_dynamics(fr, phi, tau, vcl, tcl, kprime, tauprime, tau0, model=EF_model)

        elif EF_model == 'DBM':
            peff,keff = self._calc_dynamics(fr, phi, tau, vcl, tcl, kprime, tauprime, tau0, model=EF_model)
            

        rho_tilde = mass_density+(phi*self.density_temp)-((self.density_temp**2)/peff)   #eq. 15
        
        gamma_tilde = (self.density_temp/peff)-1     #e. 16
        
        rho_eq_limp = (1/(phi*peff)) + ((gamma_tilde**2)/(phi*rho_tilde)) + (((1-phi)/phi)*(gamma_tilde/rho_tilde))  #eq. 23
        
        rho_eq_limp = 1/rho_eq_limp
        
        Zp = np.sqrt(rho_eq_limp*keff)
        kp = self.ang_freq*np.sqrt(rho_eq_limp/keff)

        TM = self._create_layer_TM(Zp,kp,thickness)
        
        return([TM,thickness,layer_name])

    
    def Add_Resistive_Screen(self,
                             thickness: float,
                             flow_resistivity: float,
                             porosity: float,
                             save_layer: bool = False,
                             layer_name: str = None) -> list[np.ndarray, float, str]:
        

        
        """
        Define a resistive screen layer

        A simplified model for thin acoustic screens -- Gaborit, Dazel, Goransson https://doi.org/10.1121/1.5047929

        Parameters
        ----------
        thickness (float):
            layer thickness [mm]
        
        flow_resisitivty (float):
            Static air flow resistivity of the material [Pa*s/m2]
            
        porosity (float):
            open porosity of the material

        save_layer (bool):
            Specify whether to save the input parameters to a database for later use.

        layer_name (str):
            If save_layer is set to True, specify the name of the layer.  Must be a unique identifier.

        Returns
        -------
        TM, thickness, layer_name (list[ndarray, float, str]):
            The transfer matrix, thickness, and name of the layer.
        
        """
        if save_layer == True:
            params = [layer_name,'screen','null',thickness,flow_resistivity,porosity,'null','null','null','null','null','null','null','null','null','null','null']
            self._layer_to_db(params)
                 
        thickness = thickness/1000
        
        fr = flow_resistivity
        phi = porosity

        w = self.ang_freq
        
        keff = self.P0/phi
        peff = ((self.density_temp/phi)+(fr/(1j*w)))

        Zp = np.sqrt(peff*keff)
        kp = self.ang_freq*np.sqrt(peff/keff)
        
        TM = self._create_layer_TM(Zp,kp,thickness)
        
        return([TM,thickness,layer_name])
    
    def Add_MAA_MPP_Layer(self,
                          thickness: float,
                          pore_diameter: float,
                          c_to_c_dist: float,
                          save_layer: bool = False,
                          layer_name: str = None) -> list[np.ndarray, float, str]:
        
        """
        Define a microperforated layer using Maa's model

        https://www.acoustics.asn.au/conference_proceedings/INTERNOISE2014/papers/p894.pdf

        Parameters
        ----------
        thickness (float):
            layer thickness [mm]
        
        pore_diameter (float):
            diameter of the microperforate [mm]
            
        c_to_c_dist (float):
            center to center distance of the microperforates [mm]

        save_layer (bool):
            Specify whether to save the input parameters to a database for later use.

        layer_name (str):
            If save_layer is set to True, specify the name of the layer.  Must be a unique identifier.

        Returns
        -------
        TM, thickness, layer_name (list[ndarray, float, str]):
            The transfer matrix, thickness, and name of the layer.
        
        """
        if save_layer == True:
            params = [layer_name,'MAA_MPP','null',thickness,'null','null','null','null','null','null','null','null','null',pore_diameter,c_to_c_dist,'null','null']
            self._layer_to_db(params)
        
        thickness = thickness/1000
        
        d = pore_diameter/1000
        b = c_to_c_dist/1000
        
        w = self.ang_freq
        
        x = d*np.sqrt((w*self.density_temp)/(4*self.viscosity_temp))
        phi = (np.pi/4)*((d/b)**2)
        
        x = d/2*np.sqrt(w*self.density_temp/(self.viscosity_temp))
        r1 = np.sqrt(1+x**2/32)+np.sqrt(2)/32*x*d/thickness
        r = 32*self.viscosity_temp/phi*thickness/d**2*r1
        m1 = 1+1/np.sqrt(1+x**2/2)+0.85*d/thickness
        m = self.density_temp*thickness/phi*m1
        
        Zp = (r+(1j*w*m))
        
        TM = self._create_Maa_MPP_TM(Zp)
        
        return([TM,thickness,layer_name])
    
    def Add_MPP_EF_Layer(self,
                         thickness: float,
                         pore_diameter: float,
                         c_to_c_dist: float,
                         save_layer: bool = False,
                         layer_name: str = None) -> list[np.ndarray, float, str]:
        """
        Define a microperforated layer using an equivalent fluid model

        Parameters
        ----------
        thickness (float):
            layer thickness [mm]
        
        pore_diameter (float):
            diameter of the microperforate [mm]
            
        c_to_c_dist (float):
            center to center distance of the microperforates [mm]
        
        save_layer (bool):
            Specify whether to save the input parameters to a database for later use.

        layer_name (str):
            If save_layer is set to True, specify the name of the layer.  Must be a unique identifier.

        Returns
        -------
        TM, thickness, layer_name (list[ndarray, float, str]):
            The transfer matrix, thickness, and name of the layer.
        
        """
        if save_layer == True:
            params = [layer_name,'EF_MPP','null',thickness,'null','null','null','null','null','null','null','null','null',pore_diameter,c_to_c_dist,'null','null']
            self._layer_to_db(params)
            
        thickness = thickness/1000
        
        d = pore_diameter/1000
        b = c_to_c_dist/1000

        phi = (np.pi/4)*((d/b)**2)

        eps = 2.0*np.sqrt(phi/np.pi)
        fok = (1-(1.13*eps)-(0.09*(eps**2))+(0.27*(eps**3)))*(4*d/3*np.pi)
        fr = (32*self.viscosity_temp)/(phi*(d**2))
        
        vcl = d/2
        tcl = d/2
        
        tau = 1+(2*fok/thickness)
        
        dyns = self._calc_dynamics(fr, phi, tau, vcl, tcl, 0, 0, 0, model='JCA')
      
        Zp = np.sqrt(dyns[0]*dyns[1])
        kp = self.ang_freq*np.sqrt(dyns[0]/dyns[1])
        
        TM = self._create_layer_TM(Zp,kp,thickness)
        
        return([TM,thickness,layer_name])

    def Add_Layer_From_Tube(self,
                            no_gap_file: str,
                            gap_file: str,
                            sample_thickness: float,
                            air_gap_thickness: float,
                            measurement: str = 'reflection') -> list[np.ndarray, float, str]:
        

        
        """
        Define a layer from the normal incidence reflection coefficients or the surface impedance of a material obtained from
        an impedance tube.  Utsuno's method currently implemented. 

        Transfer function method for measuring characteristic impedance and propagation constant of porous materials
        The Journal of the Acoustical Society of America 86, 637 (1989); https://doi.org/10.1121/1.398241

        Parameters
        ----------
        no_gap_file (str):
            Name of the csv filepath that contains the frequency dependent absorption, reflection, or surface impedance coefficients of a single porous layer obtained 
            from an impedence tube measurement with rigid backing.  The csv file should contain 2 columns of equal length -- the frequencies in the 1st column and 
            coefficients in the 2nd.
            
        gap_file (str):
            Name of the csv filepath that contains the frequency dependent absorption, reflection, or surface impedance coefficients of a single porous layer obtained 
            from an impedence tube measurement with an air gap backing.  The csv file should contain 2 columns of equal length -- the frequencies in the 1st column and 
            coefficients in the 2nd.
            
        sample_thickness (float):
            thickness of the sample [mm]

        air_gap_thickness (float):
            thickness of the air gap [mm] in the gap mounting condition.

        measurement (str):
            'reflection' or 'surface' measurement types used in the No_Gap and Gap parameters.

        Returns
        -------
        TM, thickness, layer_name (list[ndarray, float, str]):
            The transfer matrix, thickness, and name of the layer.
        
        """
        thickness = sample_thickness/1000
        air_thickness = air_gap_thickness/1000
        
        
        no_gap_data = self.load_to_array(no_gap_file)
        gap_data = self.load_to_array(gap_file)
        
        if measurement == 'surface':
            Zs_NG = no_gap_data[:,1]
            Zs_G = gap_data[:,1]

        elif measurement == 'reflection':
            Zs_NG = self.Z0*((1+no_gap_data[:,1])/(1-no_gap_data[:,1]))
            Zs_G = self.Z0*((1+gap_data[:,1])/(1-gap_data[:,1]))
        
        if np.array_equal(no_gap_data[:,0],gap_data[:,0]) != True:
            raise ValueError("Frequencies must match between no gap and gap curves")
            
        T11A = np.cos(self.k0*air_thickness)
        T21A = (1j/self.Z0)*np.sin(self.k0*air_thickness)
        Zs_A = T11A/T21A

        Zp = np.sqrt((Zs_G*(Zs_NG+Zs_A))-(Zs_NG*Zs_A))
        kp = np.arctan(Zp/(1j*Zs_NG))/thickness

        TM = self._create_layer_TM(Zp,kp,thickness)
        
        for n in range(1,10):
            rhs = int((n*self.soundspeed_temp)/(2*air_thickness))
            if rhs < self.fmax:
                print(f"Caution: results near {rhs} Hz may be unreliable due to air gap selection!")
        
        return([TM,thickness,None]) 
        
    def Add_Layer_From_Database(self,
                                layer_name: str) -> list[np.ndarray, float, str]:

        '''
        Define a layer from properties that have been saved to a database.

        Parameters
        ----------
        layer_name (str):
            The unique name the layer was saved to the database as.

        Returns
        -------
        TM, thickness, layer_name (list[ndarray, float, str]):
            The transfer matrix, thickness, and name of the layer.
        '''
        
        s = AcoustiBase()
        
        data = s.query("SELECT * from LAYER WHERE name = ?", (layer_name,))
        model_type = data[0][2]

        if model_type == 'DB':
            layer = self.Add_DB_Layer(data[0][4], data[0][5])
            return(layer)
        
        elif model_type == 'DBM':
            layer = self.Add_DBM_Layer(data[0][4], data[0][5])
            return(layer)
        
        elif model_type == 'JCA':
            layer = self.Add_JCA_Layer(data[0][4], data[0][5], data[0][6], data[0][7], data[0][8], data[0][9])
            return(layer)
        
        elif model_type == 'JCAL':
            layer = self.Add_JCAL_Layer(data[0][4], data[0][5], data[0][6], data[0][7], data[0][8], data[0][9], data[0][10])
            return(layer)
        
        elif model_type == 'JCAPL':
            layer = self.Add_JCAPL_Layer(data[0][4], data[0][5], data[0][6], data[0][7], data[0][8], data[0][9], data[0][10], data[0][11],data[0][12])
            return(layer)
        
        elif model_type == 'horoshenkov':
            layer = self.Add_Horoshenkov_Layer(data[0][4],data[0][6],data[0][16],data[0][17])
            return(layer)
        
        elif model_type == 'biot_rigid':
            layer = self.Add_Biot_Rigid_Layer(data[0][3],data[0][4], data[0][5], data[0][13], data[0][6], data[0][7], data[0][8], data[0][9], data[0][10], data[0][11], data[0][12])
            return(layer)
        
        elif model_type == 'biot_limp':
            layer = self.Add_Biot_Limp_Layer(data[0][3],data[0][4], data[0][5], data[0][13], data[0][6], data[0][7], data[0][8], data[0][9], data[0][10], data[0][11], data[0][12])
            return(layer)
        
        elif model_type == 'screen':
            layer = self.Add_Resistive_Screen(data[0][4], data[0][5], data[0][6])
            return(layer)
        
        elif model_type == 'MAA_MPP':
            layer = self.Add_MAA_MPP_Layer(data[0][4], data[0][14], data[0][15])
            return(layer)
        
        elif model_type == 'EF_MPP':
            layer = self.Add_MPP_EF_Layer(data[0][4], data[0][14], data[0][15])
            return(layer)
        
        elif model_type == 'identified':
            layer = self.Add_JCA_Layer(data[0][4], data[0][5], data[0][6], data[0][7], data[0][8], data[0][9])
            return(layer)

        elif model_type == 'AIR':
            layer = self.Add_Air_Layer(data[0][4])
            return(layer)

    def assemble_from_database(self,
                               name: str) -> list[np.ndarray, float]:
        '''
        Define a multilayer structure that has been saved to a database.

        Parameters
        ----------
        name (str):
            The unique name the multilayer structure was saved to the database as.

        Returns
        -------
        Tt,thickness (list[ndarray, float]):
            The total transfer matrix and total thickness of the structure.
        '''
        s = AcoustiBase()
        
        data = s.query("SELECT * from STRUCTURE WHERE structure_name = ?", (name,))[0]
        data = data[2:]

        layers = []

        for d in data:
            if d != 'null':
                layers.append(self.Add_Layer_From_Database(d))

        structure = self.assemble_structure(layers,db_flag=True)

        return(structure)

    def assemble_structure(self,
                           *kwargs,
                           save_structure: bool=False,
                           structure_name: str=None,
                           db_flag: bool=False) -> list[np.ndarray, float]:
        """
        Calculates the total transfer matrix for a structure of 'n' number of layers.  The structure is defined from left to right --> left being the face
        of the structure where sound impinges on the surface and right being the back or bottom of the structure that sound propagates through.

        Parameters
        ----------
        *kwargs (list):
            individual transfer matrices, thicknesses, and names of each layer, which is returned by any of the "Add_XXX_Layer" methods.

        save_structure (bool):
            Specify whether to save the structure to a database for later use.

        structure_name (str):
            If save_structure is set to True, specify the name of the structure.  Must be a unique identifier.

        db_flag (bool):
            DO NOT CHANGE THIS PARAMETER -- for internal calclations only.

        Returns
        -------
        Tt,thickness list([ndarray, float]):
            The total transfer matrix and total thickness of the structure.
        
        """

        transfer_matrices = []
        layer_thickness = []
        layer_names = []
        
        if db_flag == True:
            for i in range(len(kwargs[0])):
                transfer_matrices.append(kwargs[0][i][0])
                layer_thickness.append(kwargs[0][i][1])
                layer_names.append(kwargs[0][i][2])
                    
            
        else:
            for kwarg in kwargs:
                transfer_matrices.append(kwarg[0])
                layer_thickness.append(kwarg[1])
                layer_names.append(kwarg[2])

        
        if save_structure == True:
            s = AcoustiBase()
            data = s.pull('STRUCTURE')
            id1 = len(data)+1
            param_base = [id1,structure_name]
            for layer in layer_names:
                param_base.append(layer)
            
            params = param_base.copy()
            for i in range(18-len(param_base)):
                params.append('null')
            s.execute(params,'STRUCTURE')
            s.commit()
            s.close()
        
        thickness = sum(layer_thickness)
        if self.incidence == 'Normal':
            
            if len(transfer_matrices) == 1:
                Tt = transfer_matrices[0]
                return([Tt,thickness])
            
            elif len(transfer_matrices) > 1:
                Tt = np.einsum('ijn,jkn->ikn', transfer_matrices[0], transfer_matrices[1])
                for i in range(len(transfer_matrices)-2):
                    Tt = np.einsum('ijn,jkn->ikn', Tt, transfer_matrices[i+2])
                    
                return([Tt,thickness])
            
            elif len(transfer_matrices) == 0:
                raise ValueError("Error: Structure Not Defined. Specify each layer in assemble_structure.")
        
        elif self.incidence == "Diffuse":
            
            if len(transfer_matrices) == 1:
                Tt = transfer_matrices[0]
                return([Tt,thickness])
            
            elif len(transfer_matrices) > 1:
                Tt = np.einsum('ijnm,jknm->iknm', transfer_matrices[0], transfer_matrices[1])
                for i in range(len(transfer_matrices)-2):
                    Tt = np.einsum('ijnm,jknm->iknm', Tt, transfer_matrices[i+2])
                    
                return([Tt,thickness])
            
            elif len(transfer_matrices) == 0:
                raise ValueError("Error: Structure Not Defined. Specify each layer in assemble_structure.")

    def reflection(self,
                   transfer_matrix: list) -> np.ndarray:
        """
        Calculates the frequency dependent reflection coefficients of the structure.

        Parameters
        ----------
        transfer_matrix (list):
            total transfer matrix and thickness of the structure, which is returned by the "assemble_structure" method.

        Returns
        -------
        curve (ndarray):
            The 2D array of frequencies and reflection coefficients
        
        """
        Tt = transfer_matrix[0]
        
        if self.incidence == 'Normal':
            Zst = Tt[0][0] / Tt[1][0]
    
            R = (Zst-self.Z0)/(Zst+self.Z0)
            
            curve = np.column_stack((self.frequency,R))
            return (curve)
        
        elif self.incidence == "Diffuse":
            
            angles = np.arange(self.angles[0],self.angles[1],self.angles[2])
            v = np.cos(np.radians(angles))
            
            Zst = Tt[0][0][:][:]/ Tt[1][0][:][:]
            r = ((Zst*v)-self.Z0)/((Zst*v)+self.Z0)
            
            thetas = np.cos(np.radians(angles))*np.sin(np.radians(angles))
            
            num = np.sum(r*thetas,1)
            denom = sum(thetas)
            
            R = num/denom

            curve = np.column_stack((self.frequency,R))
            return (curve)
    
    def absorption(self,
                   transfer_matrix: list) -> np.ndarray:
        """
        Calculates the frequency dependent absorption coefficients of the structure.

        Parameters
        ----------
        transfer_matrix (list):
            total transfer matrix and thickness of the structure, which is returned by the "assemble_structure" method.

        Returns
        -------
        curve (ndarray):
            The 2D array of frequencies and absorption coefficients
        
        """
        Tt = transfer_matrix[0]
        
        if self.incidence == 'Normal':
            Zst = Tt[0][0] / Tt[1][0]
    
            R = (Zst-self.Z0)/(Zst+self.Z0)
    
            A = 1-abs(R)**2
            curve = np.column_stack((self.frequency,A))
            
            return (curve)
        
        elif self.incidence == "Diffuse":
            
            angles = np.arange(self.angles[0],self.angles[1],self.angles[2])
            v = np.cos(np.radians(angles))
            
            Zst = Tt[0][0][:][:]/ Tt[1][0][:][:]
            R = ((Zst*v)-self.Z0)/((Zst*v)+self.Z0)
            
            a = 1-abs(R)**2
            
            thetas = np.cos(np.radians(angles))*np.sin(np.radians(angles))
            
            num = np.sum(a*thetas,1)
            denom = sum(thetas)
            
            A = num/denom
            curve = np.column_stack((self.frequency,A))
            
            return (curve)
        
    def transmission_loss(self,
                          transfer_matrix: list) -> np.ndarray:
        """
        Calculates the frequency dependent transmission coefficients of the structure.

        Parameters
        ----------
        transfer_matrix (list):
            total transfer matrix and thickness of the structure, which is returned by the "assemble_structure" method.

        Returns
        -------
        curve (ndarray):
            The 2D array of frequencies and transmission coefficients
        
        """
        Tt = transfer_matrix[0]
        thickness = transfer_matrix[1]
        
        if self.incidence == 'Normal':
            
            T = (2.0*np.exp(1j*self.k0*thickness))/(Tt[0][0]+(Tt[0][1]/self.Z0)+(self.Z0*Tt[1][0])+Tt[1][1])
            Te = abs(T)**2
            TL = 10*np.log10(1/Te)

            
            curve = np.column_stack((self.frequency,TL))
            return (curve)
        
        elif self.incidence == "Diffuse":
            
            angles = np.arange(self.angles[0],self.angles[1],self.angles[2])
            v = np.cos(np.radians(angles))
            
            t1 = 2.0*np.exp(1j*self.k0*thickness)
            t1 = np.tile(t1,(int((self.angles[1]-self.angles[0])/self.angles[2]),1)).T
            t2 = (Tt[0][0][:][:]+(Tt[0][1][:][:]*v/self.Z0)+(self.Z0*Tt[1][0][:][:]/v)+Tt[1][1][:][:])
            T = t1/t2 
            Te = abs(T)**2
            
            
            thetas = np.cos(np.radians(angles))*np.sin(np.radians(angles))
            
            num = np.sum(Te*thetas,1)
            denom = sum(thetas)
            
            Tz = num/denom
            TL = 10*np.log10(1/Tz)

            curve = np.column_stack((self.frequency,TL))
            
            return (curve)
    
    def octave_bands(self,
                     curve: np.ndarray,
                     kind: str='THIRD_OCTAVE') -> np.ndarray:
        """
        Calculates the third octave or octave band absorption or transmission spectrums

        Parameters
        ----------
        curve (ndarray):
            The 2D array of frequencies and absorption or transmission coefficients, which is returend by the 'absorption' or 'transmission_loss' methods.
        
        kind (str):
            'OCTAVE' or 'THIRD_OCTAVE'
        
        Returns
        -------
        octaves (ndarray):
            The 2D array of octave bands and absorption or transmission coefficients
        
        """
        if kind == 'OCTAVE':
            
            fl = self.frequency[0]
            fu = self.frequency[-1]
            
            mid2 = 1000
            spec_mid = []  
            for i in range(int(len(self.OCTAVE_PREFERRED)/2)+1):
                mid2 = mid2/(2**(1))
            
            spec_mid = [mid2] 
            for i in range(int(len(self.OCTAVE_PREFERRED))-1):
                mid2 = mid2*(2**(1))
                spec_mid.append(mid2)
            
            spec_mid = np.asarray(spec_mid)

            bands = np.column_stack((self.OCTAVE_PREFERRED,spec_mid))
            
            lower = []
            upper = []
            for n in bands[:,1]:
                lower.append(n/((2**(1/2))**(1/1)))
                upper.append(n*((2**(1/2))**(1/1)))
            
            lower = np.asarray(lower)
            upper = np.asarray(upper)
            
            bands = np.column_stack((bands,lower))
            bands = np.column_stack((bands,upper))
            
            bands = bands[np.where((bands[:,2] > fl) & (bands[:,2] < fu))]

            octave_abs = []
            for i in range(len(bands)):    
                curve_range = curve[np.where((curve[:,0]>=bands[i,2]) & (curve[:,0]<=bands[i,3]))]
                ave = np.mean(curve_range[:,1])
                octave_abs.append(ave)
            
            octaves = np.column_stack((bands[:,0],octave_abs))
            octaves = octaves[~np.isnan(octaves).any(axis=1)]

            return(octaves)
            
        elif kind == 'THIRD_OCTAVE':
            
            fl = self.frequency[0]
            fu = self.frequency[-1]

            mid2 = 1000
            spec_mid = []  
            for i in range(int(len(self.THIRD_OCTAVE_PREFERRED)/2)+2):
                mid2 = mid2/(2**(1/3))
            
            spec_mid = [mid2] 
            for i in range(int(len(self.THIRD_OCTAVE_PREFERRED))-1):
                mid2 = mid2*(2**(1/3))
                spec_mid.append(mid2)
            
            spec_mid = np.asarray(spec_mid)

            bands = np.column_stack((self.THIRD_OCTAVE_PREFERRED,spec_mid))
            
            lower = []
            upper = []
            for n in bands[:,1]:
                lower.append(n/((2**(1/2))**(1/3)))
                upper.append(n*((2**(1/2))**(1/3)))
            
            lower = np.asarray(lower)
            upper = np.asarray(upper)
            
            bands = np.column_stack((bands,lower))
            bands = np.column_stack((bands,upper))
            
            bands = bands[np.where((bands[:,2] > fl) & (bands[:,2] < fu))]

            octave_abs = []
            for i in range(len(bands)):    
                curve_range = curve[np.where((curve[:,0]>=bands[i,2]) & (curve[:,0]<=bands[i,3]))]
                ave = np.mean(curve_range[:,1])
                octave_abs.append(ave)
            
            octaves = np.column_stack((bands[:,0],octave_abs))
            octaves = octaves[~np.isnan(octaves).any(axis=1)]
            
            return(octaves)
        
    def SAA(self,
            third_octave_curve: np.ndarray) -> float:
        """
        Calculates the average sound absorption coefficient between the 200Hz and 2500Hz third octave frequency bands.  

        Parameters
        ----------
        third_octave_curve (ndarray):
            The 2D array of frequencies and absorption, which is returend by the 'octave_bands' method.
        
        Returns
        -------
        saa (float):
            The average sound absorption coefficient, rounded to 3 decimal places
        
        """
        try:
            l = int(np.where((third_octave_curve[:,0] == 200))[0].item())
            u = int(np.where((third_octave_curve[:,0] == 2500))[0].item()+1)
            
            saa = third_octave_curve[l:u,:]
            saa = round(np.mean(saa[:,1]),3)
        except TypeError:
            raise ValueError('Unable to Calculate SAA with given frequency range!')
            return
        
        return(saa)
    
    def FFA(self,
            third_octave_curve: np.ndarray) -> float:
        """
        Calculates the four frequency average sound absorption coefficient at the 250Hz, 500Hz, 1000Hz, and 2000Hz third octave frequency bands.  

        Parameters
        ----------
        third_octave_curve (ndarray):
            The 2D array of frequencies and absorption, which is returned by the 'octave_bands' method.
        
        Returns
        -------
        ffa (float):
            The four frequency average sound absorption coefficient, rounded to 3 decimal places
        
        """
        absfreq = []
        for i in (250,500,1000,2000):
            try:
                position = int(np.where((third_octave_curve[:,0] == i))[0].item())
                absorption = third_octave_curve[position,1]
                absfreq.append(absorption)
            except TypeError:
                raise ValueError('Unable to Calculate FFA with given frequency range!')
                return
            
        ffa = round(sum(absfreq)/4,3)
        
        return(ffa)
        
    def plot_curve(self,
                   curves: list,
                   labels: list = None,
                   kind: str ='LINEAR') -> None:
        """
        Plots the frequency dependent reflection, absorption, or transmission coefficients of 1 or more structures.

        Parameters
        ----------
        curves (list):
            List of 2D arrays of frequencies and reflection, absorption, or transmission coefficients, which is retured by the 'reflection',
            'absorption', or 'transmission_loss' methods.

        labels (list):, optional
            List of strings to create a legend on the plot with labels

        kind (str):
            'LINEAR' or 'LOG' --> define whether the frequencies should be converted to a log scale for plotting purposes.
        
        """
        f, ax = plt.subplots(1)
        
        if kind == 'LINEAR':
            i=0
            for curve in curves:
                if labels is not None:
                    ax.plot(curve[:,0],curve[:,1],label=labels[i])
                    ax.legend(loc="lower right")
                    i+=1
                else:
                   ax.plot(curve[:,0],curve[:,1]) 
                   
            ax.set_ylim(bottom=0)
            plt.show()
            
            
        elif kind == 'LOG':
            i=0
            for curve in curves:
                if labels is not None:
                    logfreq = np.log10(curve[:,0])
                    ax.plot(logfreq,curve[:,1],label=labels[i])
                    ax.legend(loc="lower right")
                    i+=1
                else:
                    logfreq = np.log10(curve[:,0])
                    ax.plot(logfreq,curve[:,1])
                    
            
            ax.set_ylim(bottom=0)
            plt.show()

        return

    def to_csv(self,
               filename: str,
               data: np.ndarray) -> None:
        """
        Saves the frequency dependent reflection, absorption, or transmission coefficients of a structure to a csv file without headers.

        Parameters
        ----------
            
        filename (str):
            Name of the csv file to save data to

        data (ndarray):
            2D array of frequencies and reflection, absorption, or transmission coefficients, which is retured by the 'reflection',
            'absorption', or 'transmission_loss' methods
        
        """
        if ".csv" in filename:
            file = filename
        else:
            file = filename+".csv"

        save_path = os.path.join(file)
        np.savetxt(save_path, data, delimiter=",")

    def load_to_array(self,
                      filename: str,
                      type: str ='complex') -> None:
        """
        Loads data from csv or excel file.

        Parameters
        ----------
        filename (str):
            Name of the file to load data from

        type (str):
            type of data being loaded -- either complex or floating point data
        
        """
        if type == 'complex':
            try:
                data = np.asarray(pd.read_csv(filename,header=None).applymap(lambda s: np.complex128(s.replace('i', 'j'))))
            except Exception:
                data = np.asarray(pd.read_excel(filename,header=None).applymap(lambda s: np.complex128(s.replace('i', 'j'))))

        elif type == 'float':
            try:
                data = np.asarray(pd.read_csv(filename,header=None))
            except Exception:
                data = np.asarray(pd.read_excel(filename,header=None))  
        
        return (data)
    
    def _layer_to_db(self,
                    params: list) -> None:
            '''
            Add a layer to the database
            
            Parameters
            ----------
            params (list):
                Attributes of the given layer

            '''
            s = AcoustiBase()
            data = s.pull('LAYER')
            id1 = len(data)+1
            params.insert(0,id1)
            s.execute(params,'LAYER')
            s.commit()
            s.close()
