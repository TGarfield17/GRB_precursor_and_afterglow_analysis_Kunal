'''
Created: 25th May 2020
Here I do:

1. Get all trials for GFU 8 year background for 'Prompt+Afterglow'.
2. Doing this for the unblinding.
3. Will use these lists to obtain the local p-values.
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import astropy
import histlite as hl
import csky as cy
import os
import pandas as pd
from csky import *


import matplotlib
matplotlib.use('Agg')

ana = cy.get_analysis(cy.selections.repo,cy.selections.GFUDataSpecs.gfu_IC86)
cy.CONF['ana'] = ana

min_data_set = ana.mjd_min
max_data_set = ana.mjd_max

#Defining rh path to the directory
trials_dir =  cy.utils.ensure_dir('/data/user/kdeoskar/big_data_files/trials_files/for_unblinding/afterglow')

cy.CONF['mp_cpus'] = 5

# Load the complete GRB list
GRB_main_list = np.genfromtxt("/home/kdeoskar/projects/GRB/data/GRBdata/GRBs_withoutheader.txt", 
                               dtype = [('GRB_name','S20'),('GRB_name_Fermi','S20'),('t_trigger','S20'),
                                        ('ra','Float64'), ('decl','float'), ('pos_error','Float64'), ('T90','Float64'), 
                                        ('T90_error','Float64'),
                                        ('T90_start','S20'), ('fluence','Float64'),('fluence_error','Float64'),
                                        ('redshift','Float64'),
                                        ('T100','Float64'),
                                        ('GBM_located','bool'),
                                        ('mjd','Float64')
                                       ])


#GFU cut is between -85 and +85 degrees
#Additional cut: Excluding 14 days from either ends
my_selection = ((GRB_main_list['pos_error']>-1)*(GRB_main_list['pos_error']<0.2)*(GRB_main_list['decl']>-85.0)*(GRB_main_list['decl']<85.0)*(GRB_main_list['mjd']>(min_data_set+14.00))*(GRB_main_list['mjd']<(max_data_set-14.00)))
GRBs_of_interest = GRB_main_list[my_selection]

# we'll save, and later load, trials from this dir
bg_dir = cy.utils.ensure_dir('{}'.format(trials_dir))

def do_background_trials(ID, raGRB, decGRB, mjdGRB, N, seed):
    
    src = cy.utils.Sources(ra=np.radians(raGRB), dec=np.radians(decGRB))
    
    mjd_grb =  mjdGRB
    inj_offset = 1.0
    inj_duration = 1 / 86400.
    search_duration = 14.00 
    
    afterglow = {
        # basics
        'ana':          ana,
        'src':          src,
        'flux':         cy.hyp.PowerLawFlux(2),
        # time-dep llh stuff
        'time':         'utf',
        'box':          True,
        'box_mode':     'post',
        'dt_max':       search_duration,
        'dt_min':       0.5 * (1/86400.0),
        'seeder':       cy.seeding.UTFSeeder(threshold = 1),
        # time-dep injector
        'sig':          'tw',
        'sig_kw':       dict(t0=mjd_grb, dt=inj_offset, box_mode='post'),
        # hold t0 fixed
        'fitter_args':  dict(t0=mjd_grb),
        # use multiprocessing
        'mp_cpus':      5,
        'prior':  None #prior =None is the key for removing the penalty term
    }

    tr = cy.conf.get_trial_runner(afterglow)    
    # run trials
    trials = tr.get_many_fits(N, seed=seed, logging=False)
    # save to disk
    dir = cy.utils.ensure_dir('{}/ID/{}'.format(bg_dir, ID))
    filename = '{}/trials__N_{:06d}_seed_{:04d}_dec_{:+03.8f}.npy'.format(dir, N, seed, decGRB)
    print('->', filename)
    # notice: trials.as_array is a numpy structured array, not a cy.utils.Arrays
    np.save(filename, trials.as_array)
    
def ndarray_to_Chi2TSD(trials):
    return cy.dists.Chi2TSD(cy.utils.Arrays(trials))

#Do the background trials
background_key =5
if background_key >1:
    for x in range(len(GRBs_of_interest)):
        do_background_trials(x, GRBs_of_interest[x]['ra'], GRBs_of_interest[x]['decl'], GRBs_of_interest[x]['mjd'], 10000, seed=0)








