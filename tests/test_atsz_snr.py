from __future__ import print_function
import numpy as np
import os,sys
import pyfisher as pf
import argparse
from classy import Class
from scipy import interpolate
import subprocess
import pandas as pd

path_to_bplike = '/Users/boris/Work/CLASS-SZ/SO-SZ/bplike/'
output_dir_bplike = path_to_bplike + 'output/'
use_ps_all_freqs = True


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
    '100x143',
    '353x353',
    '217x353',
    '143x353',
    '545x545',
    '353x545',
    '217x545',
    '143x545'
             ]
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
    # list of varying foreground parameters:
    if use_ps_all_freqs:
        varying_params_fg = ['a_tsz','a_c','xi','a_ksz'] + ps_params
    else:
        varying_params_fg = ['a_tsz','a_c','a_d','a_p_tt_15','xi','a_ksz']

    varying_params = varying_params_cosmo + varying_params_fg
    print('paramaters varied in the analysis')
    print(varying_params)










    # the fiducial values of the cosmological parameters
    fiducial_cosmo_param_values = {
        '100*theta_s' : 1.042143,
        'omega_b' : 0.022032,
        'omega_cdm' : 0.12038,
        'A_s' : 2.215e-9,
        'n_s' : 0.9619,
        'tau_reio' : 0.0925
        }
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
    # sky fraction
    f_sky = 0.50

    param_list = list(fiducial_cosmo_param_values.keys())
    fg_keys = list(fiducial_foreground_param_values.keys())[:6]
    # extend param list
    param_list.extend(fg_keys)
    print('parameters used in the analysis')
    print(param_list)

    # conversion factor for the foregrounds:
    factor_fg = (2.725*1e6)**-2

    # the labels of the spectra we are considering

    frequency_list = []
    for spec in specs:
        freq1 = spec.split('x')[0]
        freq2 = spec.split('x')[0]
        if freq1 not in frequency_list:
            frequency_list.append(freq1)
        if freq2 != freq1 and freq2 not in frequency_list:
            frequency_list.append(freq2)
    frequency_list.sort()
    print('frequencies used in the analysis: ')
    print(frequency_list)


    # frequency_list = ['100','143','217','353','545']#,'857']

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
                       '100*theta_s' : fiducial_cosmo_param_values['100*theta_s'],
                       'omega_b': fiducial_cosmo_param_values['omega_b'],
                       'omega_cdm':fiducial_cosmo_param_values['omega_cdm'],
                       'A_s':fiducial_cosmo_param_values['A_s'],
                        'n_s':fiducial_cosmo_param_values['n_s'],
                       'tau_reio':fiducial_cosmo_param_values['tau_reio'],
                       # Take fixed value for primordial Helium (instead of automatic BBN adjustment)
                       'YHe':0.246,
                       }
    M = Class()
    M.set(common_settings)
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
    os.chdir(path_to_bplike)
    subprocess.call(['python',path_to_bplike+'power_atsz_analysis.py',
                     '--output_dir',output_dir_bplike,
                     '--specs',str(specs),
                     '--bplike_params_dict',str(fiducial_foreground_param_values)])
    print('bplike spectra written in : ', output_dir_bplike)



    ##########################################################
    print('computing noise power spectra')
    # compute the noise spectra
    Noise_spectra = compute_noise()
    N = Noise_spectra[1]
    # evaluate noise at effective multipoles
    for id_key,key in enumerate(list(N.keys())):
        fn = interpolate.interp1d(Noise_spectra[0],N[key])
        N[key] = fn(l_eff)

    #print(Noise_spectra[1]['217'])




    ##########################################################
    # compute covariance matrix
    print('computing data covariance matrix')
    S = np.zeros((len(frequency_list),len(frequency_list),len(l_eff)))
    for spec in specs:
        freq1 = int(spec.split('x')[0])
        freq2 = int(spec.split('x')[1])
        inu1 = frequency_list.index(str(freq1))
        inu2 = frequency_list.index(str(freq2))
        fg_dict = {}
        fg_dict['ells'], fg_dict['tot'], fg_dict['a_tsz'], fg_dict['a_c'], fg_dict['a_d'], fg_dict['a_p_tt_15'], fg_dict['xi'], fg_dict['a_ksz'] = np.loadtxt(output_dir_bplike+"spectra_l_dltt_tot_tsz_cibc_cibp_rs_tszxcib_ksz_"+str(freq1)+"_"+str(freq2)+".txt",unpack=True)



        # fg_ell = fg_D[:,0]
        binned_foregrounds = np.zeros((len(l_eff)))
        if not use_ps_all_freqs:
            for key in varying_params_fg:
                f_cl_tot = interpolate.interp1d(fg_dict['ells'],factor_fg*fg_dict[key])
                binned_f_cl_tot = f_cl_tot(l_eff)
                binned_foregrounds += binned_f_cl_tot/(l_eff*(l_eff+1.)/2./np.pi)
        else:
            for key in varying_params_fg:
                if 'ps' in key:
                    if spec in key:
                        # print('adding cl_fg for:', key)
                        f_cl_tot = interpolate.interp1d(fg_dict['ells'],factor_fg*fg_dict['a_d']+factor_fg*fg_dict['a_p_tt_15'])
                        binned_f_cl_tot = f_cl_tot(l_eff)
                        binned_foregrounds += binned_f_cl_tot/(l_eff*(l_eff+1.)/2./np.pi)
                else:
                    f_cl_tot = interpolate.interp1d(fg_dict['ells'],factor_fg*fg_dict[key])
                    binned_f_cl_tot = f_cl_tot(l_eff)
                    binned_foregrounds += binned_f_cl_tot/(l_eff*(l_eff+1.)/2./np.pi)



        S[inu1][inu2] = binned_cl_tt + binned_foregrounds
        if inu1 != inu2: S[inu2][inu1] = S[inu1][inu2]

    cov = np.zeros((len(l_eff),len(specs),len(specs)))
    # print(np.shape(S))
    # exit(0)



    for a in range(len(specs)):
        for b in range(a,len(specs)):
            speca = specs[a]
            specb = specs[b]
            freqa1 = int(speca.split('x')[0])
            freqa2 = int(speca.split('x')[1])
            freqb1 = int(specb.split('x')[0])
            freqb2 = int(specb.split('x')[1])
            inua1 = frequency_list.index(str(freqa1))
            inua2 = frequency_list.index(str(freqa2))
            inub1 = frequency_list.index(str(freqb1))
            inub2 = frequency_list.index(str(freqb2))
            i = inua1
            j = inua2
            k = inub1
            l = inub2
            cov_ixj_kxl = S[i][k]*S[j][l] + S[i][l]*S[j][k] \
            + S[i][k]*delta_iip(j,l)*N[frequency_list[j]] \
            + S[j][l]*delta_iip(i,k)*N[frequency_list[i]] \
            + S[i][l]*delta_iip(j,k)*N[frequency_list[k]] \
            + S[j][k]*delta_iip(i,l)*N[frequency_list[l]] \
            + delta_iip(i,k)*delta_iip(j,l)*N[frequency_list[i]]*N[frequency_list[j]] \
            + delta_iip(i,l)*delta_iip(j,k)*N[frequency_list[k]]*N[frequency_list[l]]
            cov_ixj_kxl  *= 1./(2.*l_eff+1.)/f_sky/delta_l
            cov[:,a,b] = cov_ixj_kxl
            if a!=b: cov[:,b,a] = cov[:,a,b].copy()




    print('computing derivatives of power spectra wrt parameters for Fisher matrix')

    # compute the derivatives with respect to comsological parameters
    # rel_tol = 5e-2
    dict_dcls = {}
    dict_binned_dcls = {}
    for key,val in fiducial_cosmo_param_values.items():
        if key == 'tau_reio':
            rel_tol = 5e-2
        else:
            rel_tol = 1e-2
        tol = val*rel_tol
        common_settings[key] = val - tol
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
        dict_dcls[key] = dcls
        f_dcl = interpolate.interp1d(cls['ell'],dcls)
        binned_dcl_tt = f_dcl(l_eff)
        dict_binned_dcls[key] = binned_dcl_tt

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
                derivs_dict[key][spec] = dict_binned_dcls[key]
                #cls_dict[key][spec] = binned_cl_tt
            else:
                freq1 = int(spec.split('x')[0])
                freq2 = int(spec.split('x')[1])
                fg_dict = {}
                fg_dict['ells'], fg_dict['tot'], fg_dict['a_tsz'], fg_dict['a_c'], fg_dict['a_d'], fg_dict['a_p_tt_15'], fg_dict['xi'], fg_dict['a_ksz'] = np.loadtxt(output_dir_bplike+"spectra_l_dltt_tot_tsz_cibc_cibp_rs_tszxcib_ksz_"+str(freq1)+"_"+str(freq2)+".txt",unpack=True)
                # fg_D = np.loadtxt(output_dir_bplike+"spectra_l_dltt_tot_tsz_cibc_cibp_rs_tszxcib_"+str(freq1)+"_"+str(freq2)+".txt")
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
        else:
            param_list_latex.append(key)

    # ['ctheta', 'ombh2', 'omch2', 'As', 'ns', 'tau',
    #                     'a_tsz','a_c','a_d','a_p_tt_15','xi','a_ksz']
    # param_list_latex = ['ctheta', 'ombh2', 'omch2', 'As', 'ns', 'tau',
    #                     'a_tsz','a_c','a_d','a_p_tt_15','xi']
    # param_list_latex = ['ctheta', 'ombh2', 'omch2', 'As', 'ns', 'tau',
    #                     'a_tsz','a_c','a_d','a_p_tt_15']
    # param_list_latex = ['ctheta', 'ombh2', 'omch2', 'As', 'ns', 'tau',
    #                     'a_tsz','a_c','a_d','a_p_tt_15']
    F = pf.FisherMatrix(Fisher,param_list_latex)
    print(F)


    sigmas = F.sigmas()
    print('marginalized errors')
    # print(sigmas)

    df = pd.DataFrame([sigmas])
    print(df)
    exit(0)

    fids = {}
    fids['ctheta'] = fiducial_cosmo_param_values['100*theta_s']
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

    for key in fids.keys():
        if key not in param_list_latex:
            fids.pop(key)




    pf.contour_plot(F,fids,'contour_atsz.pdf',name=None)





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
    Nfreqs_Planck=9 #9
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
        PS_noise_Planck[i] = 2.*(noise_arr_Planck[i] * (1.0/60.0) * (np.pi/180.0))**2.0 * np.exp( ELL*(ELL+1)* sigma_arr_Planck[i]**2. )/(2.725e6)**2
        # handle overflow due to high noise at high ell, if needed
        PS_noise_Planck[i][(np.where(PS_noise_Planck[i] > MAX_NOISE))[0]] = MAX_NOISE

    dict_results = dict(zip(freqs_Planck, PS_noise_Planck))

    return ELL,dict_results


def delta_iip(i,ip):
    if i==ip:
        result = 1.
    else:
        result = 0.
    return result



if __name__ == "__main__":
    main()
