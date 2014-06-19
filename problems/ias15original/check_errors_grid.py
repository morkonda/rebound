#!/usr/bin/env python

import subprocess
import os
import numpy as np
import matplotlib.pyplot as plt

ERRORS_FNAME = "errors_grid.txt"
N_MASS = 11
N_SIZE = 11
N_ORBITS = 10.
MASS_MIN = 1.e-3
MASS_MAX = 1.e3
SIZE_MIN = 1.e-3
SIZE_MAX = 1.e3
REL_EN_ERR = dict(ndx= 7,text="Rel. Energy Error",text2="en_err")
REL_KE_ERR = dict(ndx= 8,text="Rel. KE Error",text2="ke_err")
REL_PE_ERR = dict(ndx= 9,text="Rel. PE Error",text2="pe_err")
REL_VV_ERR = dict(ndx=11,text="Rel. Vel Error",text2="v_err")
IAS15_eps = 100.
EVERHART_eps = 1.e-4

def call_rebound_grid(mass_scale_vec,size_scale_vec,
                      ias15new=True,rel_errors_list=[REL_EN_ERR]):
    errors_A = {}
    errors_filename = {}
    version = "ias15"
    if ias15new:
        rebound_error_fname = "energy_ias15variable.txt"
    else:
        rebound_error_fname = "energy_ias15original.txt"
        version += "original"
    for rel_error in rel_errors_list:
        key = rel_error["text2"]
        this_A = np.zeros((len(mass_scale_vec),
                           len(size_scale_vec)))
        errors_A[key] = np.NaN * this_A
        this_err_fname = ERRORS_FNAME[:-4] + '_' + version \
                         + '_' + key + ERRORS_FNAME[-4:]
        errors_filename[key] = this_err_fname
        try: # remove this error file if it exists
            os.remove(this_err_fname)
        except OSError: # swallow error and do nothing
            pass
    # compile (make) binary
    make_command = "make %s"%version
    subprocess.call(make_command,shell=True)
    for mass_scale_ndx in xrange(len(mass_scale_vec)):
        mass_scale = mass_scale_vec[mass_scale_ndx]
        for size_scale_ndx in xrange(len(size_scale_vec)):
            size_scale = size_scale_vec[size_scale_ndx]
            try: # remove the REBOUND error file if it exists
                os.remove(rebound_error_fname)
            except OSError: # swallow error and do nothing
                pass
            runstring = "./nbody"
            runstring += " --scale_m=%8.4e"%mass_scale
            runstring += " --scale_a=%8.4e"%size_scale
            runstring += " --tmax_factor=%4.2f"%N_ORBITS
            if ias15new: runstring += " --epsilon=%4.2e"%IAS15_eps
            else:        runstring += " --epsilon=%4.2e"%EVERHART_eps
            print runstring
            subprocess.call(runstring,shell=True)
            # harvest errors from the rebound run
            rebound_errors = np.loadtxt(rebound_error_fname)
            for rel_error in rel_errors_list:
                key = rel_error["text2"]
                this_error = rebound_errors[rel_error["ndx"]]
                errors_A[key][mass_scale_ndx,size_scale_ndx] = this_error
                with open(errors_filename[key],'a') as errors_f:
                    errors_f.write("%8.4e %8.4e %8.4e\n"%(mass_scale,
                                                          size_scale,
                                                          this_error))
    return dict(errors_A=errors_A,rel_errors_list=rel_errors_list)

def generate_grid(ias15new=True,rel_errors_list=[REL_EN_ERR]):
    mass_scale_vec = np.logspace(np.log10(MASS_MIN),np.log10(MASS_MAX),N_MASS)
    size_scale_vec = np.logspace(np.log10(SIZE_MIN),np.log10(SIZE_MAX),N_SIZE)
    sizes_A,masses_A = np.meshgrid(size_scale_vec,mass_scale_vec)
    grid_dict = call_rebound_grid(mass_scale_vec,size_scale_vec,
                                  ias15new=ias15new,rel_errors_list=rel_errors_list)
    return dict(masses_A = masses_A,
                sizes_A = sizes_A,
                grid_dict = grid_dict)

def display_grid(info_dict,ias15new=True):
    masses_A = info_dict["masses_A"]
    mass_vec = masses_A[:,0]
    sizes_A  = info_dict["sizes_A"]
    size_vec = sizes_A[0,:]
    # set up grid-cell edges vectors
    logmass_vec_edges = np.zeros(len(mass_vec)+1)
    logsize_vec_edges = np.zeros(len(mass_vec)+1)
    logmass_vec_edges[:-1] = np.log10(mass_vec) - 0.5*(np.log10(mass_vec[1])
                                                       - np.log10(mass_vec[0]))
    logmass_vec_edges[-1] = np.log10(mass_vec[-1]) + 0.5*(np.log10(mass_vec[1])
                                                          - np.log10(mass_vec[0]))
    mass_vec_edges = 10.**logmass_vec_edges
    logsize_vec_edges[:-1] = np.log10(size_vec) - 0.5*(np.log10(size_vec[1])
                                                       - np.log10(size_vec[0]))
    logsize_vec_edges[-1] = np.log10(size_vec[-1]) + 0.5*(np.log10(size_vec[1])
                                                          - np.log10(size_vec[0]))
    size_vec_edges = 10.**logsize_vec_edges
    # loop over the different errors
    grid_dict = info_dict["grid_dict"]
    errors_A = grid_dict["errors_A"]
    rel_errors_list = grid_dict["rel_errors_list"]
    for rel_error in rel_errors_list:
        key = rel_error["text2"]
        this_err = errors_A[key]
        log_abs_err = np.log10(np.abs(this_err) + 1.00001e-17)
        plt.ion()
        plt.pcolor(size_vec_edges,mass_vec_edges,log_abs_err)
        plt.xlabel("Size Scale",size=12)
        plt.ylabel("Mass Scale",size=12)
        plt.gca().set_xscale('log')
        plt.gca().set_yscale('log')
        plt.colorbar()
        plt.autoscale(tight=True)
        plt.clim(-11,-4)
        title_str = "log$_{10}$[|" + rel_error["text"] + "|] "
        if ias15new:
            title_str += "(IAS15 time-stepping)"
            plt.title(title_str,size=16)
            plt.show()
            plot_fname = "errors_new_" + key + ".pdf"
            plt.savefig(plot_fname)
            plt.close()
        else:
            title_str += "(Everhart time-stepping)"
            plt.title(title_str,size=16)
            plt.show()
            plot_fname = "errors_original_" + key + ".pdf"
            plt.savefig(plot_fname)
            plt.close()

if __name__=="__main__":
    ias15new_list = [False,True]
    rel_errors_list = [REL_EN_ERR,REL_KE_ERR,REL_PE_ERR]
    for ias15new_bool in ias15new_list:
        ret_dict = generate_grid(ias15new=ias15new_bool,
                                 rel_errors_list=rel_errors_list)
        display_grid(ret_dict,ias15new=ias15new_bool)




# old

# def read_grid(ias15new=True):
#     version = "ias15"
#     if ias15new: pass
#     else: version += "original"
#     errors_filename = ERRORS_FNAME[:-4] + '_' + version + ERRORS_FNAME[-4:]
#     data = np.loadtxt(errors_filename)
#     masses = data[:,0]
#     sizes  = data[:,1]
#     errors = data[:,2]
#     new_shape = (N_MASS,N_SIZE)
#     return dict(masses_A = np.reshape(masses,new_shape),
#                 sizes_A  = np.reshape(sizes,new_shape),
#                 errors_A = np.reshape(errors,new_shape))

        # make_grid = True
        # #make_grid = False
        # if make_grid:
        #     ret_dict = generate_grid(ias15new=ias15new_bool,rel_errors_list)
        # else: # read grid from file
        #     ret_dict = read_grid(ias15new=ias15new_bool)
