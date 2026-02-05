"""
acoustipy - Tools for characterizing the acoustic performance of porous materials.

This package provides implementations of the acoustic transfer matrix method and
optimization routines for identifying material parameters from impedance tube measurements.
"""

import torch

from acoustipy.TMM import AcousticTMM
from acoustipy.Params import AcousticID
from acoustipy.Database import AcoustiBase

__all__ = ["AcousticTMM", "AcousticID", "AcoustiBase", "cuda_available", "get_device"]


def cuda_available() -> bool:
    """
    Check if CUDA is available for GPU acceleration.
    
    Returns
    -------
    bool
        True if CUDA is available, False otherwise.
    """
    return torch.cuda.is_available()


def get_device(prefer_cuda: bool = True) -> str:
    """
    Get the best available device for computation.
    
    Parameters
    ----------
    prefer_cuda : bool, optional
        If True, return 'cuda' when available. Default is True.
        
    Returns
    -------
    str
        'cuda' if CUDA is available and preferred, otherwise 'cpu'.
    """
    if prefer_cuda and torch.cuda.is_available():
        return 'cuda'
    return 'cpu'
