import numpy as np
import matplotlib.pyplot as plt
import scipy as sp



def plot_style():

    # Set the style preferences:
    plt.style.use("seaborn-v0_8-darkgrid")
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 14,
        "axes.labelsize": 16,
        "axes.titlesize": 18,
        "legend.fontsize": 12,
        "figure.figsize": (8,6),
        "axes.grid": True,
        "grid.alpha": 0.3,
        "grid.linestyle": "--",
        "text.usetex": True,
        "xtick.labelsize": 12,      
        "ytick.labelsize": 12,      
        "xtick.major.size": 6,      
        "xtick.minor.size": 3,      
        "ytick.major.size": 6,      
        "ytick.minor.size": 3,      
        "xtick.major.width": 1.0,   
        "xtick.minor.width": 0.75,  
        "ytick.major.width": 1.0,   
        "ytick.minor.width": 0.75,  
        "xtick.direction": "out",   
        "ytick.direction": "out",   
        "xtick.color": "black",     
        "ytick.color": "black",
    })

