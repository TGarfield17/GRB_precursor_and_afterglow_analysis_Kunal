'''
Date created: 30th April 2020

This file is to get the BP distribution. We give 'Ntrials' and a 'save' key which is used to save as well as used as the seed.

This includes the dt_min change to 0.5 s as recommended by Ignacio.
'''

import argparse, os, time
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import astropy
import histlite as hl
import csky as cy
import os
import pandas as pd
import random 
from astropy.coordinates import SkyCoord
import pickle

from csky.ipyconfig import *
from csky import *
import matplotlib
matplotlib.use('Agg')

p = argparse.ArgumentParser(description="Getting BP distribution for scrambled datasets.",formatter_class=argparse.RawTextHelpFormatter)
p.add_argument("--N_trials", default=100, type=int,
                help="Number of scrambled trials")
p.add_argument("--savekey", default=1.0, type=float,
                help="Name for save key for the data files")
args = p.parse_args()

# set up timer
timer = cy.timing.Timer()
time = timer.time

def BinomialTest(pValues, kmax, returnArray=False):
    pSorted = np.sort(pValues)
    n = len(pSorted)
    if (kmax==0):
        kmax = n       # default:  search the whole array for best p-value
    i = np.arange(0,n)
    # NB the array index i = k-1.  We just need i in next line since we want sf(k-1) = sf(k)+pmf(k)
    BpArr = stats.binom.sf(i[:kmax],n,pSorted[:kmax])  # remember: test is w.r.t. n, even if we only search up to kmax
    BpBest = BpArr.min()
    kBest = BpArr.argmin()+1
    pThresh = pSorted[kBest-1]
    if returnArray:
        return pThresh, kBest, BpBest, BpArr
    else:
        return pThresh, kBest, BpBest

    
def RandomizedTrials(n, kmax, trials):
    pThreshArr = np.ones(trials)
    kBestArr = np.ones(trials)
    BpBestArr = np.ones(trials)
    for i in np.arange(trials):
        pThreshArr[i], kBestArr[i], BpBestArr[i] = BinomialTest(np.random.random(n),kmax)
    return pThreshArr, kBestArr, BpBestArr

def ndarray_to_Chi2TSD(trials):
    return cy.dists.Chi2TSD(cy.utils.Arrays(trials))

repo = selections.Repository()
with time('ana setup'):
    ana = cy.get_analysis(cy.selections.repo,cy.selections.GFUDataSpecs.gfu_IC86)
    cy.CONF['ana'] = ana
    
min_data_set = ana.mjd_min
max_data_set = ana.mjd_max
# Load the complete GRB list
GRB_main_list = np.genfromtxt("/home/kdeoskar/projects/GRB/data/GRBdata/GRBs_withoutheader.txt", 
                               dtype = [('GRB_name','S20'),('GRB_name_Fermi','S20'),('t_trigger','S20'),
                                        ('ra','f4'), ('decl','float'), ('pos_error','f4'), ('T90','f4'), ('T90_error','f4'),
                                        ('T90_start','S20'), ('fluence','f4'),('fluence_error','f4'),
                                        ('redshift','f4'),
                                        ('T100','f4'),
                                        ('GBM_located','bool'),
                                        ('mjd','f4')
                                       ])

#GFU cut is between -85 and +85 degrees
my_selection = ((GRB_main_list['pos_error']>-1)*(GRB_main_list['pos_error']<0.2)*
                (GRB_main_list['decl']>-85.0)*(GRB_main_list['decl']<85.0)*
                (GRB_main_list['mjd']>(min_data_set+14.00))*(GRB_main_list['mjd']<(max_data_set-14.00)))
GRBs_of_interest = GRB_main_list[my_selection]

trials_dir = cy.utils.ensure_dir('/data/user/kdeoskar/big_data_files/trials_files/trials_GFU/2011_2018/10k_each')
bg_dir = cy.utils.ensure_dir('{}/bg'.format(trials_dir))

with time('load bg trials'):
    bg = cy.bk.get_all(
        # disk location
        '{}/ID'.format(bg_dir),
        # filename pattern
        'trials*npy',
        # how to combine items within each directory
        merge=np.concatenate,
        # what to do with items after merge
        post_convert=ndarray_to_Chi2TSD)

TS_above_zero_list = []
for ID in range(0,len(bg)):
    count = 0
    for ts_temp in bg[ID].trials['ts']:
        if ts_temp > 0.0 : 
            count = count + 1
    TS_above_zero_list.append(count)
    
TS_above_zero_fraction = np.array(TS_above_zero_list)/len(bg[0].trials['ts'])

sources = cy.sources(GRBs_of_interest['ra'], GRBs_of_interest['decl'], deg=True)
MJD_GRB = GRBs_of_interest['mjd'] #+(GRBs_of_interest['T100']/86400.0)

#Defining the tr
def give_me_tr(raGRB, decGRB, mjd_GRB):
    
    src = cy.utils.Sources(ra=np.radians(raGRB), dec=np.radians(decGRB))
    
    mjd_grb =  mjd_GRB
    #The next two variables are only for the injection cases
    inj_offset = 0.5 
    inj_duration = 30 / 86400.
    
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
        'sig_kw':       dict(t0=mjd_grb + inj_offset, dt=inj_duration, box_mode='pre'),
        #Note that the box_mode here is set to 'pre'. This is so that we can control injecting the size of pulses.
        # hold t0 fixed
        'fitter_args':  dict(t0=mjd_grb),
        # use multiprocessing
        'mp_cpus':      10,
        'prior':  None #prior =None is the key for removing the penalty term
    }
    #print(raGRB)
    return cy.conf.get_trial_runner(afterglow, extra_keep=['run', 'event'], cut_n_sigma=5)#Made the change here on 08/04/2020

with time('trs setup'):
    trs = [give_me_tr(grb['ra'], grb['decl'], grb['mjd']
                      #+(grb['T100']/86400.0)
                     ) for grb in GRBs_of_interest]

#defining the injector

#inj_offset = 0.5 
inj_duration = 14.0
mjd_grb = 58386.61

afterglow = {
        # basics
        'ana':          ana,
        'src':          sources,
        'flux':         cy.hyp.PowerLawFlux(2),
        # time-dep llh stuff
        'time':         'utf',
        'box':          True,
        'box_mode':     'post',
        'dt_max':       14.0,
        'dt_min':       0.5 * (1/86400.0),
        'seeder':       cy.seeding.UTFSeeder(threshold = 1),
        # time-dep injector
        'inj_conf':     dict(src=sources),
        'sig':          'tw',
        'sig_kw':       dict(t0=mjd_grb, dt=inj_duration, box_mode='post'),
        #Note that the box_mode here is set to 'pre'. This is so that we can control injecting the size of pulses.
        # hold t0 fixed
        'fitter_args':  dict(t0=mjd_grb),
        # use multiprocessing
        'mp_cpus':      10,
        'prior':  None #prior =None is the key for removing the penalty term
        }

afterglow_inj = {}
afterglow_inj.update(afterglow)
del afterglow_inj['fitter_args']

tr_inj = cy.get_trial_runner(afterglow_inj, src=sources, inj_conf=dict(src=sources),extra_keep=['run', 'event']
                            )

#Defining the Multi trial runner to perform the analysis on the same scrambled sky
multr = cy.trial.MultiTrialRunner(
    # the Analysis
    ana,
    # bg+sig injection trial runner (produces trials)
    tr_inj,
    # llh test trial runners (perform fits given trials)
    trs,
    # background distrubutions, can also do without this info since was getting some error when i use it
    #bgs=bg,
    # use multiprocessing
    mp_cpus=10,
)

#Now we do the real thingy
N_TRIALS = args.N_trials
SAVE = args.savekey

print('Starting the analysis now')
print('This is for'+str(N_TRIALS)+' trials and '+str(SAVE)+' key')

with time('Many fits Multi tr'):
    M_trials = multr.get_many_fits(N_TRIALS, n_sig=0, seed=int(SAVE))

Temp_MTR = M_trials.as_dataframe
Temp_MTR.to_pickle("/home/kdeoskar/projects/GRB/analysis_csky/GRB_pickle_data/complete_analysis/BP_distribution/"+str(N_TRIALS)+"_"+str(SAVE)+"_afterglow.pkl")

unpickled_Temp_MTR = pd.read_pickle("/home/kdeoskar/projects/GRB/analysis_csky/GRB_pickle_data/complete_analysis/BP_distribution/"+str(N_TRIALS)+"_"+str(SAVE)+"_afterglow.pkl")

GRB_Multi_TS = np.zeros((733,N_TRIALS))
GRB_Multi_T_w = np.zeros((733,N_TRIALS))

#To get the data file from where we will get the new TS fraction

with time('Converting dataframes'):
    for x in range(len(GRB_Multi_TS)):
        temp_str = str('ts_')+str('{:04d}'.format(x))
        GRB_Multi_TS[x] = unpickled_Temp_MTR[temp_str]
        temp_str_tw = str('dt_')+str('{:04d}'.format(x))
        GRB_Multi_T_w[x] = unpickled_Temp_MTR[temp_str_tw]


GRB_Multi_PV = np.zeros((733,N_TRIALS))
Ts_local_Mtr_1 = GRB_Multi_TS[:,0]

with time('TS to p-v'):   
    for Y in range(len(GRB_Multi_PV[0])): #to go columnwise
        Ts_local_Mtr_1 = GRB_Multi_TS[:,Y]
        p_v_local = np.zeros(len(GRBs_of_interest))
        
        for X in range(len(GRBs_of_interest)):
                ts_local = Ts_local_Mtr_1[X]
                if ts_local >0.0:
                    background_TS = bg[X]
                    p_value_of_interest = background_TS.sf(ts_local)
                else:
                    p_value_of_interest = 1.0
                p_v_local[X] = p_value_of_interest
        GRB_Multi_PV[:,Y] = p_v_local

list_of_k_best = np.zeros(len(GRB_Multi_PV[0]))
list_of_p_local_best = np.zeros(len(GRB_Multi_PV[0]))
list_of_BP_best = np.zeros(len(GRB_Multi_PV[0]))
list_of_p_local_best = np.zeros(len(GRB_Multi_PV[0]))

with time('BT'): 
    for Y in range(len(GRB_Multi_PV[0])):#to go columnwise
        
        pValues_here = GRB_Multi_PV[:,Y]
        kmax_localmax = len(pValues_here)
        pBest, kBest, BpBest, BpArr = BinomialTest(pValues_here,kmax = kmax_localmax, returnArray=True)
        list_of_BP_best[Y] = BpBest
        list_of_k_best[Y] = kBest
        list_of_p_local_best[Y] = pBest
        
filename_BPB = "/home/kdeoskar/projects/GRB/analysis_csky/GRB_pickle_data/complete_analysis/BP_distribution/"+str(N_TRIALS)+"_"+str(SAVE)+"_BPP_list.p"
pickle.dump(list_of_BP_best, open(filename_BPB, "wb" ) )
    