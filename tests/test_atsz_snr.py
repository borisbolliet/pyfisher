from __future__ import print_function
import numpy as np
import os,sys
import pyfisher as pf
import argparse
from classy import Class
from scipy import interpolate
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import math

# COMMAND:
# $ python tests/test_atsz_snr.py

path_to_bplike = '/Users/boris/Work/CLASS-SZ/SO-SZ/bplike/'
output_dir_bplike = path_to_bplike + 'output/'
use_ps_all_freqs = False # known/unknown frequency dependence of Poisson components
add_tau_prior = 'yes'
tau_prior_sigma = 7e-3 # Planck 2018 constraint
do_contours = 'no'
plot_spectra_and_noise = 'yes'
FIG_DIR = '/Users/boris/Work/CLASS-SZ/SO-SZ/figures'
plot_covmat =  False


# Here is what this file does:
# Computation of primary cmb with class
# Computation of foreground spectra with bplike (Mat M. implementation)
# Computation of noise power spectra (jch implementation for Planck)
# Computation of data covariance matrix
# Computation of Fisher matrix (based on pyfisher implementation from Mat M.)
# Contour plot for free-varying parameters


def parse_args(args):
    # Parse command line
    parser = argparse.ArgumentParser(description='Do a thing.')
    parser.add_argument("--exp", type=str,default='planck',help='Experiment name.')
    parser.add_argument("--Lmin",     type=int,  default=4,help="A description.")
    parser.add_argument("--Lmax",     type=int,  default=400,help="A description.")
    parser.add_argument("--fsky",     type=float,  default=0.65,help="A description.")
    args = parser.parse_args()
    return args

def main():
    args = parse_args(sys.argv[1:])
    run_fisher_forecast(args.Lmin,args.Lmax,args.exp,args.fsky,test=False)

def run_fisher_forecast(Lmin,Lmax,exp,fsky,test=False):


    print('starting computation')
    specs = [
    '100x100',
    '143x143',
    '217x217',
    '143x217',
    # '353x217'

    # '353x353',
    # '217x353',
    # '143x353',
    #
    # '545x545',
    # '353x545',
    # '217x545',
    # '143x545',
    #
    # '100x143',
    # '100x217',
    # '100x353',
    # '100x545',
             ]

    # specs = [
    #         'TT',
    #         'EE',
    #         'TE'
    #         ]


    # specs = ['100x100', '143x143', '217x217', '143x217', '100x143', '353x353', '217x353', '143x353', '545x545', '353x545', '217x545', '143x545']

    print('considering the following components:')
    print(specs)
    spec_type_TEB = bool([True for s in specs if 'T' in s or 'E' in s or 'B' in s])
    if spec_type_TEB:
        print('spec type TEB, not currently setup to use fg')
    # exit(0)
    else:
        if use_ps_all_freqs:
            print('definning params')
            ps_params = []
            for spec in specs:
                pstr = 'a_ps_' + spec
                ps_params.append(pstr)
        # print(ps_params)
    # exit(0)





    # list of varying cosmological parameters:
    varying_params_cosmo = ['100*theta_s','omega_b','omega_cdm','A_s','n_s','tau_reio']
    # varying_params_cosmo = ['H0','omega_b','omega_cdm','A_s','n_s','tau_reio']
    # list of varying foreground parameters:
    if spec_type_TEB:
        varying_params_fg = []
    else:
        if use_ps_all_freqs:
            varying_params_fg = ['a_tsz','a_c','xi','a_ksz'] + ps_params
        else:
            varying_params_fg = ['a_tsz','a_c','a_d','a_p_tt_15' ,'xi','a_ksz']
        # varying_params_fg = []


    varying_params = varying_params_cosmo + varying_params_fg
    print('paramaters varied in the analysis')
    print(varying_params)










    # the fiducial values of the cosmological parameters
    fiducial_cosmo_param_values = {
        '100*theta_s' : 1.042143,
        # 'H0' : 67.02393,
        'omega_b' : 0.022032,
        'omega_cdm' : 0.12038,
        'A_s' : 2.215e-9,
        'n_s' : 0.9619,
        'tau_reio' : 0.0925
        }
    # fiducial_cosmo_param_values = {
    #     'H0' : 67.02393,
    #     'omega_b' : 0.02219218,
    #     'omega_cdm' : 0.1203058,
    #     'A_s' : 2.15086031154146e-9,
    #     'n_s' : 0.9625356,
    #     'tau_reio' : 0.06574325
    #     }

    # sky fraction
    f_sky = 0.60
    # f_sky = 1.

    # param_list = list(fiducial_cosmo_param_values.keys())
    # fg_keys = list(fiducial_foreground_param_values.keys())[:6]
    # extend param list
    # print(fg_keys)
    # param_list.extend(fg_keys)
    # print('parameters used in the analysis')
    # print(param_list)
    comp_dict = {}
    comp_list = []
    # conversion factor for the foregrounds:
    factor_fg = (2.725*1e6)**-2
    # if not spec_type_TEB:
    #
    #
    #     # the labels of the spectra we are considering
    #
    #     comp_list = []
    #     for spec in specs:
    #         comp1 = spec.split('x')[0]
    #         comp2 = spec.split('x')[1]
    #         print(spec,comp1,comp2)
    #         if comp1 not in comp_list:
    #             comp_list.append(comp1)
    #         if comp2 not in comp_list:
    #             comp_list.append(comp2)
    #     comp_list.sort()
    #     print('frequencies used in the analysis: ')
    #     print(comp_list)
    # else:
    for spec in specs:
        if spec_type_TEB:
            comp1,comp2 = spec
        else:
            comp1 = spec.split('x')[0]
            comp2 = spec.split('x')[1]
        print(comp1,comp2)
        if comp1 not in comp_list:
            comp_list.append(comp1)
        if comp2 not in comp_list:
            comp_list.append(comp2)
    comp_list.sort()
    print('components used in the analysis: ')
    print(comp_list)
    for id_comp,comp in enumerate(comp_list):
        comp_dict[id_comp] = comp
    print(comp_dict)
    # exit(0)


    # exit(0)


    # comp_list = ['100','143','217','353','545']#,'857']

    # exit(0)




    ##########################################################

    # define the binning
    l_min = 10 # minimum ell
    l_max = 2508 # maximum ell
    delta_l = 20 # number of ell's per bin

    nbin = int(np.floor((l_max - l_min)/delta_l))

    l_low = []
    l_p = l_min
    ib = 0
    while ib < nbin:
        l_low.append(l_p)
        l_p += delta_l
        ib += 1
    l_low = np.asarray(l_low)
    l_up = l_low[1:] - 1
    l_up = np.append(l_up,l_low[-1]+delta_l-1)

    # effective multipoles:
    l_eff = (l_low + l_up) / 2

    print('binning used in the analysis: l_min = %d, l_max = %d, dl = %d, nbin = %d'%(l_min,l_max,delta_l,nbin))

    # factor_binned = 1.e10*l_eff*(l_eff+1.)/2./math.pi

    ##########################################################
    # compute primary cl's fiducial using class:
    common_settings = {# wich output? ClTT, transfer functions delta_i and theta_i
                       'output':'tCl',
                       # LambdaCDM parameters
                       #'h':0.67556,
                       # 'H0' : fiducial_cosmo_param_values['H0'],

                       # '100*theta_s' : fiducial_cosmo_param_values['100*theta_s'],
                       # 'omega_b': fiducial_cosmo_param_values['omega_b'],
                       # 'omega_cdm':fiducial_cosmo_param_values['omega_cdm'],
                       # 'A_s':fiducial_cosmo_param_values['A_s'],
                       #  'n_s':fiducial_cosmo_param_values['n_s'],
                       # 'tau_reio':fiducial_cosmo_param_values['tau_reio'],
                       # Take fixed value for primordial Helium (instead of automatic BBN adjustment)
                       'YHe':0.246,
                       }
    for key in list(fiducial_cosmo_param_values.keys()):
        common_settings[key] = fiducial_cosmo_param_values[key]

    M = Class()
    M.set(common_settings)

    if not spec_type_TEB:
        M.set({'output':'tCl','modes':'s','lensing':'no','l_max_scalars':5000})
        M.compute()
        cls = M.raw_cl(5000)
        M.struct_cleanup()
        M.empty()
        # binned cl_tt
        f_cl = interpolate.interp1d(cls['ell'],cls['tt'])
        binned_cl_tt = f_cl(l_eff)
        # print(binned_cl_tt)

        # print(str(fiducial_foreground_param_values))
        # exit(0)

        ##########################################################
        # compute foreground cl's fiducial using bplike:
        print('Running bplike')
        # the fiducial values of the foreground parameters
        fiducial_foreground_param_values = {
            # first the foreground amplitudes in this order:
            'a_tsz': 3.3,
            'a_c': 4.9,
            'a_d': 6.9,
            'a_p_tt_15': 3.1,
            'xi': 1.654936e-2,
            'a_ksz':  4.950210,
            # then other bplike paramaters
            'flux_cut': '15mJy'
            }
        os.chdir(path_to_bplike)
        subprocess.call(['python',path_to_bplike+'power_atsz_analysis.py',
                         '--output_dir',output_dir_bplike,
                         '--specs',str(specs),
                         '--bplike_params_dict',str(fiducial_foreground_param_values)])
        print('bplike spectra written in : ', output_dir_bplike)
    else:
        pol = bool([True for s in specs if 'E' in s or 'B' in s])
        if pol:
            M.set({'output':'tCl,pCl','modes':'s','lensing':'no','l_max_scalars':5000})
        else:
            M.set({'output':'tCl,pCl','modes':'s','lensing':'no','l_max_scalars':5000})
        M.compute()
        cls = M.raw_cl(5000)
        M.struct_cleanup()
        M.empty()
        binned_cl_tt = {}
        for id_comp1 in range(len(comp_list)):
            for id_comp2 in range(id_comp1,len(comp_list)):
                comp1 = comp_dict[id_comp1]
                comp2 = comp_dict[id_comp2]
                s = comp1+comp2
                s_key = s.lower()
                try:
                    f_cl = interpolate.interp1d(cls['ell'],cls[s_key])
                except KeyError:
                    f_cl = interpolate.interp1d(cls['ell'],cls[s_key[::-1]])
                # f_cl = interpolate.interp1d(cls['ell'],cls['tt'])
                binned_cl_tt[s] = f_cl(l_eff)
                binned_cl_tt[s[::-1]] = binned_cl_tt[s].copy()


        # for s in specs:
        #     s_key = s.lower()
        #     # f_cl = interpolate.interp1d(cls['ell'],cls[s_key])
        #     f_cl = interpolate.interp1d(cls['ell'],cls['tt'])
        #     binned_cl_tt[s] = f_cl(l_eff)
        #     binned_cl_tt[s[::-1]] = binned_cl_tt[s].copy()
    # print(binned_cl_tt)
    # exit(0)


    ##########################################################
    print('computing noise power spectra')
    nls_dict = {}
    if not spec_type_TEB:
        # compute the noise spectra
        Noise_spectra = compute_noise()
        # Noise_spectra = compute_noise_mat()
        N_dict = Noise_spectra[1]
        # evaluate noise at effective multipoles
        N = {}
        # print('enum',list(N.keys()))
        for id_key,key in enumerate(list(N_dict.keys())):
            fn = interpolate.interp1d(Noise_spectra[0],N_dict[key])
            N[key] = fn(l_eff)
            # print('noise:',key,N[key])


        for id_keya,keya in enumerate(comp_list):
            nls_dict_key = str(id_keya)+str(id_keya)
            nls_dict[nls_dict_key] = N[keya]
    else:
        Noise_spectra = pf.get_planck_nls(l_eff)
        for id_keya,keya in comp_dict.items():
            nls_dict_key = str(id_keya)+str(id_keya)
            nls_dict[nls_dict_key] = Noise_spectra[keya+keya](l_eff)/(2.725e6)**2
            # nls_dict[nls_dict_key] = Noise_spectra['TT'](l_eff)/(2.725e6)**2


    # print(nls_dict)
    # exit(0)

        #print(Noise_spectra[1]['217'])





    ##########################################################
    # compute covariance matrix
    print('computing data covariance matrix')
    S = np.zeros((len(comp_list),len(comp_list),len(l_eff)))

    cls_dict = {}
    binned_fg_dict = {}
    for id_comp1 in range(len(comp_list)):
        for id_comp2 in range(id_comp1,len(comp_list)):
            if not spec_type_TEB:
                spec = comp_list[id_comp1] + 'x' + comp_list[id_comp2]
                binned_fg_dict[spec] = {}
                comp1 = int(spec.split('x')[0])
                comp2 = int(spec.split('x')[1])

                inu1 = comp_list.index(str(comp1))
                inu2 = comp_list.index(str(comp2))


                fg_dict = {}
                fg_dict['ells'], fg_dict['tot'], fg_dict['a_tsz'], fg_dict['a_c'], fg_dict['a_d'], fg_dict['a_p_tt_15'], fg_dict['xi'], fg_dict['a_ksz'] = np.loadtxt(output_dir_bplike+"spectra_l_dltt_tot_tsz_cibc_cibp_rs_tszxcib_ksz_"+str(comp1)+"_"+str(comp2)+".txt",unpack=True)



                # fg_ell = fg_D[:,0]
                binned_foregrounds = np.zeros((len(l_eff)))
                if not use_ps_all_freqs:
                    for key in varying_params_fg:
                        f_cl_tot = interpolate.interp1d(fg_dict['ells'],factor_fg*fg_dict[key])
                        binned_f_cl_tot = f_cl_tot(l_eff)
                        binned_foregrounds += binned_f_cl_tot/(l_eff*(l_eff+1.)/2./np.pi)
                        binned_fg_dict[spec][key] = binned_f_cl_tot/(l_eff*(l_eff+1.)/2./np.pi)
                else:
                    for key in varying_params_fg:
                        if 'ps' in key:
                            if spec in key:
                                # print('adding cl_fg for:', key)
                                f_cl_tot = interpolate.interp1d(fg_dict['ells'],factor_fg*fg_dict['a_d']+factor_fg*fg_dict['a_p_tt_15'])
                                binned_f_cl_tot = f_cl_tot(l_eff)
                                binned_foregrounds += binned_f_cl_tot/(l_eff*(l_eff+1.)/2./np.pi)
                                binned_fg_dict[spec][key] = binned_f_cl_tot/(l_eff*(l_eff+1.)/2./np.pi)
                        else:
                            f_cl_tot = interpolate.interp1d(fg_dict['ells'],factor_fg*fg_dict[key])
                            binned_f_cl_tot = f_cl_tot(l_eff)
                            binned_foregrounds += binned_f_cl_tot/(l_eff*(l_eff+1.)/2./np.pi)
                            binned_fg_dict[spec][key] = binned_f_cl_tot/(l_eff*(l_eff+1.)/2./np.pi)


                S[inu1,inu2,:] = binned_cl_tt + binned_foregrounds
            else:
                # comp1,comp2 = spec
                inu1 = id_comp1 #comp_list.index(str(comp1))
                inu2 = id_comp2 #comp_list.index(str(comp2))
                spec = comp_dict[inu1]+comp_dict[inu2]
                S[inu1,inu2,:] = binned_cl_tt[spec]

            if inu1 != inu2: S[inu2,inu1,:] = S[inu1,inu2,:].copy()
            cls_key = str(inu1) + str(inu2)
            cls_dict[cls_key] = S[inu1,inu2,:].copy()
            cls_dict[cls_key[::-1]] = cls_dict[cls_key].copy()

    # print(cls_dict)
    # exit(0)




    cov = np.zeros((len(l_eff),len(specs),len(specs)))
    # print(np.shape(S))
    # exit(0)


    #
    # for a in range(len(specs)):
    #     for b in range(a,len(specs)):
    #         speca = specs[a]
    #         specb = specs[b]
    #         freqa1 = int(speca.split('x')[0])
    #         freqa2 = int(speca.split('x')[1])
    #         freqb1 = int(specb.split('x')[0])
    #         freqb2 = int(specb.split('x')[1])
    #         inua1 = comp_list.index(str(freqa1))
    #         inua2 = comp_list.index(str(freqa2))
    #         inub1 = comp_list.index(str(freqb1))
    #         inub2 = comp_list.index(str(freqb2))
    #         i = inua1
    #         j = inua2
    #         k = inub1
    #         l = inub2
    #         cov_ixj_kxl = S[i][k]*S[j][l] + S[i][l]*S[j][k] \
    #         + S[i][k]*delta_iip(j,l)*N[comp_list[j]] \
    #         + S[j][l]*delta_iip(i,k)*N[comp_list[i]] \
    #         + S[i][l]*delta_iip(j,k)*N[comp_list[k]] \
    #         + S[j][k]*delta_iip(i,l)*N[comp_list[l]] \
    #         + delta_iip(i,k)*delta_iip(j,l)*N[comp_list[i]]*N[comp_list[j]] \
    #         + delta_iip(i,l)*delta_iip(j,k)*N[comp_list[k]]*N[comp_list[l]]
    #         cov_ixj_kxl  *= 1./(2.*l_eff+1.)/f_sky/delta_l
    #         cov[:,a,b] = cov_ixj_kxl
    #         if a!=b:
    #             print('symmetrizing')
    #             cov[:,b,a] = cov[:,a,b].copy()

    # for a in range(len(specs)):
    #     for b in range(a,len(specs)):
    #         speca = specs[a]
    #         specb = specs[b]
    #         print(speca,specb)
    #         freqa1 = int(speca.split('x')[0])
    #         freqa2 = int(speca.split('x')[1])
    #         freqb1 = int(specb.split('x')[0])
    #         freqb2 = int(specb.split('x')[1])
    #         print(freqa1,freqa2,freqb1,freqb2)
    #         inua1 = comp_list.index(str(freqa1))
    #         inua2 = comp_list.index(str(freqa2))
    #         inub1 = comp_list.index(str(freqb1))
    #         inub2 = comp_list.index(str(freqb2))
    #         i = inua1
    #         j = inua2
    #         k = inub1
    #         l = inub2
    #
    #
    #
    #         cl_ik = S[i,k]
    #         # print(a,b,i,k,delta_iip(i,k))
    #         nl_ik =  delta_iip(i,k)*N[comp_list[i]]
    #         cl_jl = S[j,l]
    #         nl_jl =  delta_iip(j,l)*N[comp_list[j]]
    #         cl_il = S[i,l]
    #         nl_il =  delta_iip(i,l)*N[comp_list[l]]
    #         cl_jk = S[j,k]
    #         nl_jk =  delta_iip(j,k)*N[comp_list[k]]
    #         # print('i,j,k,l :',i,j,k,l)
    #         # print('n_ik : ', nl_ik)
    #         # print('n_jl : ', nl_jl)
    #         # print('n_il : ', nl_il)
    #         # print('n_jk : ', nl_jk)
    #
    #
    #         cov[:,a,b] = ((cl_ik+nl_ik)*(cl_jl+nl_jl)+(cl_il+nl_il)*(cl_jk+nl_jk))/(2*l_eff+1)/delta_l/f_sky
    #         if a!=b: cov[:,b,a] = cov[:,a,b].copy()

    # Mat's version
    def _symmz(cdict,ab):
        try:
            return cdict[ab]
        except KeyError:
            try:
                return cdict[ab[::-1]]
            except KeyError:
                # print('problem in nls/cls_dict with key: ', ab)
                return np.zeros((len(l_eff),))


    ncomps = len(specs)
    nbins = len(l_eff)
    cov = np.zeros((nbins,ncomps,ncomps))
    for i in range(ncomps):
        for j in range(i,ncomps):
            spec1 = specs[i]
            spec2 = specs[j]
            # print(spec1,spec2)
            # print(type(spec1))
            if not spec_type_TEB:
                spec1_key = str(comp_list.index(str(int(spec1.split('x')[0]))))+str(comp_list.index(str(int(spec1.split('x')[1]))))
                spec2_key = str(comp_list.index(str(int(spec2.split('x')[0]))))+str(comp_list.index(str(int(spec2.split('x')[1]))))
            else:
                spec1_key = str(comp_list.index(str(spec1[0])))+str(comp_list.index(str(spec1[1])))
                spec2_key = str(comp_list.index(str(spec2[0])))+str(comp_list.index(str(spec2[1])))

            a,b = spec1_key
            g,d = spec2_key



            ag = a+g
            bd = b+d
            ad = a+d
            bg = b+g


            cl_ag = _symmz(cls_dict,ag)
            nl_ag = _symmz(nls_dict,ag)
            cl_bd = _symmz(cls_dict,bd)
            nl_bd = _symmz(nls_dict,bd)
            cl_ad = _symmz(cls_dict,ad)
            nl_ad = _symmz(nls_dict,ad)
            cl_bg = _symmz(cls_dict,bg)
            nl_bg = _symmz(nls_dict,bg)
            # print('a,b :',a,b)
            # print('g,d :',g,d)
            # print('a,g,nl_ag :', a,g,nl_ag*(2.725e6)**2)
            # print('a,g,cl_ag :', a,g,cl_ag*(2.725e6)**2)
            # print('b,d,nl_bd :', b,d,nl_bd*(2.725e6)**2)
            # print('b,d,cl_bd :', b,d,cl_bd*(2.725e6)**2)
            # print('a,d,nl_ad :', a,d,nl_ad*(2.725e6)**2)
            # print('a,d,cl_ad :', a,d,cl_ad*(2.725e6)**2)
            # print('b,g,nl_bg :', b,g,nl_bg*(2.725e6)**2)
            # print('b,g,cl_bg :', b,g,cl_bg*(2.725e6)**2)

            cov[:,i,j] = ((cl_ag+nl_ag)*(cl_bd+nl_bd)+(cl_ad+nl_ad)*(cl_bg+nl_bg))/(2*l_eff+1)/delta_l/f_sky
            if i!=j: cov[:,j,i] = cov[:,i,j].copy()


    # exit(0)
    # print('condition covmat')
    # print(np.linalg.cond(cov))
    # print('cov[-1,:,:]')
    # print(cov[-1,:,:] )
    print('np.shape(cov): ', np.shape(cov))
    print(specs)

    if plot_covmat:
        # covmat for plotting:

        covmat = np.zeros((len(l_eff)*len(specs),len(l_eff)*len(specs)))
        print(covmat)
        for a in range(len(specs)):
            for b in range(len(specs)):
                for ell in range(len(l_eff)):
                    covmat[(a-1)*len(l_eff)+ell][(b-1)*len(l_eff)+ell] = cov[ell,a,b]

        print(covmat)
        # def extents(f):
        #   delta = f[1] - f[0]
        #   return [f[0] - delta/2, f[-1] + delta/2]
        plt.imshow(np.log10(covmat.transpose()), aspect='auto', interpolation='none',origin='upper')
        plt.savefig(FIG_DIR +'py.png')

    # cobmat = np.block([[A11, A12, .. , A1N], [A21, A22, .., A2N]])
    # print(np.shape(cov))
    # print('cov[0,:,:]')
    # print(cov[0,:,:] )


    print('computing derivatives of power spectra wrt parameters for Fisher matrix')

    # compute the derivatives with respect to comsological parameters
    # rel_tol = 5e-2
    # dict_dcls = {}
    dict_binned_dcls = {}
    for key,val in fiducial_cosmo_param_values.items():
        if key == 'tau_reio':
            rel_tol = 5e-2
        else:
            rel_tol = 1e-3
        tol = val*rel_tol
        common_settings[key] = val - tol
        if not spec_type_TEB:
            M = Class()
            M.set(common_settings)
            M.set({'output':'tCl','modes':'s','lensing':'no','l_max_scalars':5000})
            M.compute()
            cls = M.raw_cl(5000)
            M.struct_cleanup()
            M.empty()
            cls_low = cls['tt']

            common_settings[key] = val + tol
            M = Class()
            M.set(common_settings)
            M.set({'output':'tCl','modes':'s','lensing':'no','l_max_scalars':5000})
            M.compute()
            cls = M.raw_cl(5000)
            M.struct_cleanup()
            M.empty()
            cls_up = cls['tt']

            dcls = (cls_up - cls_low)/2./tol
            # dict_dcls[key] = dcls
            f_dcl = interpolate.interp1d(cls['ell'],dcls)
            binned_dcl_tt = f_dcl(l_eff)
            dict_binned_dcls[key] = binned_dcl_tt
        else:
            pol = bool([True for s in specs if 'E' in s or 'B' in s])
            M = Class()
            M.set(common_settings)
            if pol:
                M.set({'output':'tCl,pCl','modes':'s','lensing':'no','l_max_scalars':5000})
            else:
                M.set({'output':'tCl','modes':'s','lensing':'no','l_max_scalars':5000})
            M.compute()
            cls = M.raw_cl(5000)
            M.struct_cleanup()
            M.empty()
            cls_low = {}
            for s in specs:
                s_key = s.lower()
                try:
                    cls_low[s] = cls[s_key]
                except KeyError:
                    cls_low[s] = cls[s_key[::-1]]
                # cls_low[s] = cls['tt']

            common_settings[key] = val + tol
            M = Class()
            M.set(common_settings)
            if pol:
                M.set({'output':'tCl,pCl','modes':'s','lensing':'no','l_max_scalars':5000})
            else:
                M.set({'output':'tCl','modes':'s','lensing':'no','l_max_scalars':5000})
            M.compute()
            cls = M.raw_cl(5000)
            M.struct_cleanup()
            M.empty()
            cls_up = {}
            dict_binned_dcls[key] = {}
            for s in specs:
                s_key = s.lower()
                cls_up[s] = cls[s_key]
                # cls_up[s] = cls['tt']

                dcls = (cls_up[s] - cls_low[s])/2./tol
                #dict_dcls[key] = dcls
                f_dcl = interpolate.interp1d(cls['ell'],dcls)
                binned_dcl_tt = f_dcl(l_eff)
                dict_binned_dcls[key][s] = binned_dcl_tt

    # print('derivatives cosmo params')
    # print(dict_binned_dcls)




    # blike notations (from act_pylike.py)
    # a_tsz : t_sz amplitude
    # a_c : clustered cib
    # a_d : poisson cib
    # xi : tsz_x_cib
    # a_p_tt_15 : # TT radio Poisson with given flux cut in mJy





    derivs_dict = {}
    #cls_dict = {}
    for key in varying_params:
        derivs_dict[key] = {}
        for spec in specs:
            if key not in varying_params_fg:
                if not spec_type_TEB:
                    derivs_dict[key][spec] = dict_binned_dcls[key]
                else:
                    derivs_dict[key][spec] = dict_binned_dcls[key][spec]
                #cls_dict[key][spec] = binned_cl_tt
            else:
                comp1 = int(spec.split('x')[0])
                comp2 = int(spec.split('x')[1])
                fg_dict = {}
                fg_dict['ells'], fg_dict['tot'], fg_dict['a_tsz'], fg_dict['a_c'], fg_dict['a_d'], fg_dict['a_p_tt_15'], fg_dict['xi'], fg_dict['a_ksz'] = np.loadtxt(output_dir_bplike+"spectra_l_dltt_tot_tsz_cibc_cibp_rs_tszxcib_ksz_"+str(comp1)+"_"+str(comp2)+".txt",unpack=True)
                # fg_D = np.loadtxt(output_dir_bplike+"spectra_l_dltt_tot_tsz_cibc_cibp_rs_tszxcib_"+str(comp1)+"_"+str(comp2)+".txt")
                if not use_ps_all_freqs:
                    f_cl_tot = interpolate.interp1d(fg_dict['ells'],factor_fg*fg_dict[key])
                    binned_f_cl_tot = f_cl_tot(l_eff)
                    derivs_dict[key][spec] = binned_f_cl_tot/(l_eff*(l_eff+1.)/2./np.pi)
                else:
                    if 'ps' not in key:
                        f_cl_tot = interpolate.interp1d(fg_dict['ells'],factor_fg*fg_dict[key])
                        binned_f_cl_tot = f_cl_tot(l_eff)
                        derivs_dict[key][spec] = binned_f_cl_tot/(l_eff*(l_eff+1.)/2./np.pi)
                    else:
                        if spec in key:
                            f_cl_tot = interpolate.interp1d(fg_dict['ells'],factor_fg*fg_dict['a_d']+factor_fg*fg_dict['a_p_tt_15'])
                            binned_f_cl_tot = f_cl_tot(l_eff)
                            derivs_dict[key][spec] = binned_f_cl_tot/(l_eff*(l_eff+1.)/2./np.pi)
                        else:
                            derivs_dict[key][spec] = 0.




    #
    # print('full dictionnar of dcls')
    # print(derivs_dict)
    print('computing Fisher matrix')
    cinv = np.linalg.inv(cov)
    nparams = len(varying_params)
    Fisher = np.zeros((nparams,nparams))

    for i in range(nparams):
        for j in range(i,nparams):

            param1 = varying_params[i]
            param2 = varying_params[j]
            dcls1 = np.zeros((len(l_eff),len(specs)))
            dcls2 = np.zeros((len(l_eff),len(specs)))
            for k,spec in enumerate(specs):
                dcls1[:,k] = derivs_dict[param1][spec]
                dcls2[:,k] = derivs_dict[param2][spec]

            Fisher[i,j] = np.einsum('ik,ik->',np.einsum('ij,ijk->ik',dcls1,cinv),dcls2)
            if i!=j: Fisher[j,i] = Fisher[i,j]




    param_list_latex = []
    for key in varying_params:
        if key == '100*theta_s':
            param_list_latex.append('ctheta')
        elif key == 'omega_b':
            param_list_latex.append('ombh2')
        elif key == 'omega_cdm':
            param_list_latex.append('omch2')
        elif key == 'A_s':
            param_list_latex.append('As')
        elif key == 'n_s':
            param_list_latex.append('ns')
        elif key == 'tau_reio':
            param_list_latex.append('tau')
        else:
            param_list_latex.append(key)


    F = pf.FisherMatrix(Fisher,param_list_latex)
    if 'tau' in param_list_latex and add_tau_prior == 'yes':
        F.add_prior('tau',tau_prior_sigma)
    print(F)


    sigmas = F.sigmas()
    print('marginalized errors')
    print('a_tsz :',sigmas['a_tsz'])
    print('a_ksz :',sigmas['a_ksz'])
    # print(sigmas)
    pd.set_option('display.float_format', lambda x: '%.2e' % x)
    df = pd.DataFrame([sigmas])
    print(df)
    if do_contours == 'yes':
        fids = {}
        try:
            fids['ctheta'] = fiducial_cosmo_param_values['100*theta_s']
        except:
            fids['H0'] = fiducial_cosmo_param_values['H0']
        fids['ombh2'] = fiducial_cosmo_param_values['omega_b']
        fids['omch2'] = fiducial_cosmo_param_values['omega_cdm']
        fids['As'] = fiducial_cosmo_param_values['A_s']
        fids['ns'] = fiducial_cosmo_param_values['n_s']
        fids['tau'] = fiducial_cosmo_param_values['tau_reio']
        fids['a_tsz'] = 1.0
        fids['a_c'] = 1.0
        fids['a_d'] = 1.0
        fids['a_p_tt_15'] = 1.0
        fids['xi'] = 1.0
        fids['a_ksz'] = 1.0

        fids_key = list(fids.keys()).copy()
        for key in fids_key:
            if key not in param_list_latex:
                fids.pop(key)

        pf.contour_plot(F,fids,'contour_atsz.pdf',name=None)
    if plot_spectra_and_noise == 'yes':
        nplots = len(specs)
        fig, ax = plt.subplots(nplots,1,figsize=(7,5*nplots))
        for ids,spec in enumerate(specs):
            try:
                ax_id = ax[ids]
            except TypeError:
                ax_id = ax
            ax_id.set_xlim([90,2000])
            #ax_id.set_ylim([-50,50])
            ax_id.set_xlabel(r"$\ell$")
            ax_id.set_ylabel(r"$\ell (\ell+1) C_l^{XY} / 2 \pi \,\,\, [\times 10^{10}]$")
            factor_binned = 1.e10*l_eff*(l_eff+1.)/2./math.pi
            if not spec_type_TEB:
                ax_id.plot(l_eff,factor_binned*binned_cl_tt,marker='o',
                           markersize = 3,label='primary',ls='-')
                comp1 = int(spec.split('x')[0])
                comp2 = int(spec.split('x')[1])
                inu1 = comp_list.index(str(comp1))
                inu2 = comp_list.index(str(comp2))
                ax_id.plot(l_eff,factor_binned*S[inu1][inu2] ,marker='o',
                           markersize = 3,label='primary+fg',ls='-')

                yerr = factor_binned*np.sqrt(cov[:,ids,ids])
                ax_id.errorbar(l_eff,factor_binned*binned_cl_tt,yerr=yerr,label='cl+cov')
                for key in varying_params_fg:
                    if 'ps' in key:
                        if spec in key:
                            ax_id.plot(l_eff,np.abs(factor_binned*binned_fg_dict[spec][key]),label=key)
                    else:
                        ax_id.plot(l_eff,np.abs(factor_binned*binned_fg_dict[spec][key]),label=key)
            else:
                yerr = factor_binned*np.sqrt(cov[:,ids,ids])
                ax_id.errorbar(l_eff,factor_binned*binned_cl_tt[spec],yerr=yerr,label='cl+cov')
                ax_id.errorbar(l_eff,-factor_binned*binned_cl_tt[spec],yerr=yerr,ls='--')

            ax_id.legend(loc='right',bbox_to_anchor=(1.4, 0.5),fontsize=10)
            ax_id.set_yscale('log')
            ax_id.set_xscale('linear')
            ax_id.set_title(spec)
            ax_id.grid(which='both')
        plt.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.25,
                            wspace=0.35)
        fig.tight_layout()
        plt.show(block=False)
        FIG_NAME = '/cl_XY'
        plt.savefig(FIG_DIR + FIG_NAME +".pdf")
        print('plot saved')


    exit(0)


    #
    # # Load default fiducial parameters
    # fids = pyfisher.get_fiducials()
    #
    # # Load a pre-calculated BOSS BAO Fisher
    # bao = pyfisher.get_saved_fisher('boss_bao')
    #
    # # Load a pre-calculated CMB lensing noise curve
    # ells,nls = pyfisher.get_lensing_nl(exp)
    # # Calculate a CMB lensing Fisher
    # bin_edges = np.arange(Lmin,Lmax)
    # lens = pyfisher.get_lensing_fisher(bin_edges,ells,nls,fsky)
    #
    # # Planck lens + BAO (s8, om, H0 parameterization)
    # F = lens+bao
    # F.delete(['w0','wa','ok','nnu','tau','mnu'])
    # s8 = pyfisher.get_s8(zs=[0.],params=fids)[0]
    # fids['s8'] = s8
    # fids['om'] = (fids['omch2'] + fids['ombh2'])/(fids['H0']/100)**2.
    # F1 = pyfisher.reparameterize(F,['om','s8','H0','ns','ombh2'],fids,verbose=False)
    # F1.add_prior('ns',0.02)
    # F1.add_prior('ombh2',0.0005)
    # sigmas = F1.sigmas()
    # print("Planck lens + BAO (s8, om, H0 parameterization)")
    # for p in ['s8','om','H0']:
    #     print(f'{p} = {fids[p]:.03f}+-{sigmas[p]:.03f}')
    # if test:
    #     assert np.isclose(sigmas['s8'],0.01806,rtol=1e-3)
    #     assert np.isclose(sigmas['om'],0.02242,rtol=1e-3)
    #
    # # Planck lens alone (s8om^0.25, H0 parameterization)
    # fids['s8om0.25'] = fids['s8'] * fids['om']**0.25
    # F = pyfisher.reparameterize(lens,['s8om0.25','H0','ns','ombh2'],fids,verbose=False)
    # F.add_prior('ns',0.02)
    # F.add_prior('ombh2',0.0005)
    # sigmas = F.sigmas()
    # print("Planck lens  (s8om^0.25, H0 parameterization)")
    # for p in ['s8om0.25']:
    #     print(f'{p} = {fids[p]:.03f}+-{sigmas[p]:.03f}')
    # if test:
    #     assert np.isclose(sigmas['s8om0.25'],0.01879,rtol=1e-3)
    #
    #
    # # Planck lens + BAO + CMB (mnu)
    # lcmb = pyfisher.get_saved_fisher('planck_lowell',0.65)
    # hcmb = pyfisher.get_saved_fisher('planck_highell',0.65)
    # F = lens + bao + lcmb + hcmb
    # F.delete(['w0','wa','ok','nnu'])
    # F.add_prior('tau',0.007)
    # sigmas = F.sigmas()
    # print("Planck lens + BAO + CMB (mnu)")
    # for p in ['mnu']:
    #     print(f'{p} = {fids[p]:.03f}+-{sigmas[p]:.03f}')
    # if test:
    #     assert np.isclose(sigmas['mnu'],0.0797,rtol=1e-3)
    # else:
    #     pyfisher.contour_plot(F,fids,'contour.png',name=None)
#
# def test_lensing_demo():
#     run_lensing_demo(4,400,'planck',0.65,test=False)

def compute_noise():
    # Planck freqs -- no 545 or 857
    Nfreqs_Planck= 9 #9
    freqs_Planck = []
    freqs_Planck.append('030')
    freqs_Planck.append('044')
    freqs_Planck.append('070')
    freqs_Planck.append('100')
    freqs_Planck.append('143')
    freqs_Planck.append('217')
    freqs_Planck.append('353')
    freqs_Planck.append('545')
    freqs_Planck.append('857')
    #freqs_Planck_float = np.array([30.0e9, 44.0e9, 70.0e9, 100.0e9, 143.0e9, 217.0e9, 353.0e9])
    freqs_Planck_float = np.array([30.0e9, 44.0e9, 70.0e9, 100.0e9, 143.0e9, 217.0e9, 353.0e9,545.0e9,857.0e9])

    # Planck noise
    noise_arr_Planck = np.zeros(Nfreqs_Planck)
    noise_arr_Planck[0] = 195.079975053 #uK-arcmin, from Table 7 (first column) of https://arxiv.org/pdf/1502.01585.pdf -- converted via sqrt(3224.4*(4*Pi*(180/Pi)^2*60^2/(12*1024^2)))
    noise_arr_Planck[1] = 226.090506617 # converted via sqrt(4331.0*(4*Pi*(180/Pi)^2*60^2/(12*1024^2)))
    noise_arr_Planck[2] = 199.09525581 # converted via sqrt(3358.5*(4*Pi*(180/Pi)^2*60^2/(12*1024^2))) assuming Nside=1024... I think this is the right assumption (would be lower for 2048)
    noise_arr_Planck[3] = 77.4 #uK-arcmin, from Table 6 of https://arxiv.org/pdf/1502.01587v2.pdf
    noise_arr_Planck[4] = 33.0
    noise_arr_Planck[5] = 46.8
    noise_arr_Planck[6] = 153.6
    noise_arr_Planck[7] = 0.78 #kJy/sr * deg -- need to convert to uK-arcmin
    noise_arr_Planck[8] = 0.72 #kJy/sr * deg -- need to convert to uK-arcmin

    # Planck beams
    FWHM_arr_Planck = np.zeros(Nfreqs_Planck)
    FWHM_arr_Planck[0] = 32.239
    FWHM_arr_Planck[1] = 27.005
    FWHM_arr_Planck[2] = 13.252
    FWHM_arr_Planck[3] = 9.69 #arcmin, from Table 6 of https://arxiv.org/pdf/1502.01587v2.pdf
    FWHM_arr_Planck[4] = 7.30
    FWHM_arr_Planck[5] = 5.02
    FWHM_arr_Planck[6] = 4.94
    FWHM_arr_Planck[7] = 4.83
    FWHM_arr_Planck[8] = 4.64
    # convert to sigma in radians
    sigma_arr_Planck = FWHM_arr_Planck / np.sqrt(8. * np.log(2)) /60. * np.pi/180.

    # Planck noise power spectra
    ELLMAX=8000
    ELL = np.arange(ELLMAX+1)
    MAX_NOISE = 1.e9 #prevent overflow due to high noise in Planck
    PS_noise_Planck = np.zeros((Nfreqs_Planck,ELLMAX+1))
    for i in range(Nfreqs_Planck):#factpr 2 for half mission
        PS_noise_Planck[i] = (noise_arr_Planck[i] * (1.0/60.0) * (np.pi/180.0))**2.0 * np.exp( ELL*(ELL+1)* sigma_arr_Planck[i]**2. )/(2.725e6)**2

        # handle overflow due to high noise at high ell, if needed
        PS_noise_Planck[i][(np.where(PS_noise_Planck[i] > MAX_NOISE))[0]] = MAX_NOISE

    dict_results = dict(zip(freqs_Planck, PS_noise_Planck))

    return ELL,dict_results


def compute_noise_mat():
    Nfreqs_Planck=7 #9
    freqs_Planck = []
    freqs_Planck.append('030')
    freqs_Planck.append('044')
    freqs_Planck.append('070')
    freqs_Planck.append('100')
    freqs_Planck.append('143')
    freqs_Planck.append('217')
    freqs_Planck.append('353')
    # Planck freqs -- no 545 or 857
    beams_T =  [33.,23.,14.,10.,7.,5.,5.]
    uK_arcmins_T = [145.,149.,137.,65.,43.,66.,200.]
    beams_P =  [14.,10.,7.,5.,5.]
    uK_arcmins_P = [450.,103.,81.,134.,406.]
    ELLMAX=3000
    ELL = np.arange(ELLMAX+1)
    ells = ELL
    Ns_TT = np.asarray([(uK_arcmin*np.pi/180./60.)**2./pf.gauss_beam(ells,fwhm)**2. for uK_arcmin,fwhm in zip(uK_arcmins_T,beams_T)])
    # print('shape of ns_TT: ',np.shape(Ns_TT))
    # Ns_PP = np.asarray([(uK_arcmin*np.pi/180./60.)**2./gauss_beam(ells,fwhm)**2. for uK_arcmin,fwhm in zip(uK_arcmins_P,beams_P)])
    # N_TT = 1./(1./Ns_TT).sum(axis=0)
    # N_PP = 1./(1./Ns_PP).sum(axis=0)
    PS_noise_Planck = np.zeros((Nfreqs_Planck,ELLMAX+1))
    for i in range(Nfreqs_Planck):#factpr 2 for half mission
        PS_noise_Planck[i] = Ns_TT[i]/(2.725e6)**2#2.*(noise_arr_Planck[i] * (1.0/60.0) * (np.pi/180.0))**2.0 * np.exp( ELL*(ELL+1)* sigma_arr_Planck[i]**2. )/(2.725e6)**2

        # handle overflow due to high noise at high ell, if needed
        # PS_noise_Planck[i][(np.where(PS_noise_Planck[i] > MAX_NOISE))[0]] = MAX_NOISE

    dict_results = dict(zip(freqs_Planck, PS_noise_Planck))

    return ELL,dict_results



def delta_iip(i,ip):
    if i==ip:
        result = 1.
    else:
        result = 0.
    return result


#
# def get_planck_nls(ells):
#     beams_T =  [33.,23.,14.,10.,7.,5.,5.]
#     uK_arcmins_T = [145.,149.,137.,65.,43.,66.,200.]
#     beams_P =  [14.,10.,7.,5.,5.]
#     uK_arcmins_P = [450.,103.,81.,134.,406.]
#     Ns_TT = np.asarray([(uK_arcmin*np.pi/180./60.)**2./gauss_beam(ells,fwhm)**2. for uK_arcmin,fwhm in zip(uK_arcmins_T,beams_T)])
#     # print('shape of ns_TT: ',np.shape(Ns_TT))
#     Ns_PP = np.asarray([(uK_arcmin*np.pi/180./60.)**2./gauss_beam(ells,fwhm)**2. for uK_arcmin,fwhm in zip(uK_arcmins_P,beams_P)])
#     N_TT = 1./(1./Ns_TT).sum(axis=0)
#     N_PP = 1./(1./Ns_PP).sum(axis=0)
#     nls = {}
#     # print('noise ell, nl TT')
#     # print(ells)
#     # print(N_TT)
#     nls['TT'] = interp(ells,N_TT)
#     nls['EE'] = interp(ells,N_PP)
#     nls['BB'] = interp(ells,N_PP)
#     return nls



if __name__ == "__main__":
    main()
