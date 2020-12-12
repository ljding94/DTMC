#!/usr/local/bin/python3
import time
import numpy as np
import matplotlib.pyplot as plt
from plot import *
from Oplot import *
from check_plot import *
from analyze import *


def main():
    print("hello!")
    # Shape_helix_energy()
    # curvature_compare()
    #config_plot_xyz("../data/Ne2/Jul25_2020/State_N200_Ne2_L0_kar10_C00.0_karg0.0_lam6.0_B10.0_Bc0_tau00.00.txt",tag=r"$\lambda=6.0,\bar{\kappa}=0.0$")
    #config_plot3D("../data/scratch_local/State_N200_Ne1_L10_kar0_lam0.0_karg0.0_B0.0_Bc0_tau00.00_init.txt")
    #config_plot3D("../data/Ne2/May17_2020/State_N200_Ne1_L0_kar10_lam2.0_karg2.0_B10.0_Bc0_tau00.00.txt")


    #config_plot_xyz("../State_error.txt")
    foldername = "../data/Ne2/Aug1_2020"
    print("analyzing "+foldername)
    pars = []
    colors = []
    alphas = []
    Ns = [100,200,300,400]
    N = 100
    Ne=2
    Nes = [1,2]
    Ls = [13,14,15, 16, 17,18, 19, 20, 21, 22, 23, 24, 25, 26,27, 28, 30, 31, 32, 33, 35, 36, 37, 38, 39, 40,41,42,43,44,45,46,47,48,49,50]
    # run with longer L range
    L = 0
    kars = [1,5,10,15]
    #kars = [0.2,0.4,0.6,0.8,1.0]
    kar = 40
    C0=0.0
    karg=0.0
    kargs=[0.0,1.0,2.0,3.0,4.0]
    #kargs=np.linspace(-3.0,3.0,13)
    #kargs = [-1.0,-0.5,0.0,0.5,1.0]
    #Bs = np.linspace(0.0,15.0,31)
    lam =4.0
    #lams = np.linspace(0.0,11.0,56)
    lams = np.linspace(2.0,8.0,13)
    #lams = [2.0]
    Bs = np.linspace(0.0,8.0,5)
    B = 10.0
    Bc=0
    Bcs=[20,40,60,80,100,200]
    tau0=0.00
    #tau0s = [-0.10,-0.20,0.00,0.10,0.20,0.30]
    #tau0s = np.linspace(0.00,0.30,7)
    #tau0s = [0.00,0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90]
    tau0s = [0.00,0.10,0.20,0.30,0.40,0.50,0.60,0.70]
    #tau0s = [0.00,0.10,0.20,0.30,0.40]
    Ne1pars = []
    Ne2pars = []
    for karg in kargs:
        #Ne1pars.append([N, 1, L, kar, lam, karg, B, Bc, tau0])
        #Ne2pars.append([N, 2, L, kar, lam, karg, B, Bc, tau0])
        pars.append([N, Ne, L, kar, C0, karg, lams, B, Bc, tau0])
    #pars.append([N, Ne, Ls, kar, lam, 4.0, B, Bc, 0.00])
    colors = None
    alphas = None
    scatterfilenames = []
    for i in range(len(pars)):
        N, Ne, L, kar, C0, karg,lam, B, Bc, tau0 = pars[i]
        print(pars[i])
        O_stat_ana(foldername, N,Ne, L, kar,C0, karg,lam, B, Bc, tau0, 6, mode="lam")
        #Cxx_stat_plot(foldername,"ncnc", N,Ne, L, kar, lam,karg, B, Bc, tau0, "lam",2)
        #Cxx_stat_plot(foldername,"rr", N,Ne, L, kar, lam,karg, B, Bc, tau0, "tau0",1)
        # -- config plot

        for lam in lams:
            pass
            filename = foldername + "/State_N%.0f_Ne%.0f_L%.0f_kar%.0f_C0%.1f_karg%.1f_lam%.1f_B%.1f_Bc%.0f_tau0%.2f.txt" % (N,Ne, L, kar, C0, karg, lam,B, Bc, tau0)
            #Fix_index = [int(N / (2 * L)) * L-1,int(N / (2 * L)) * L-1+L-1]
            if int(karg)==karg and int(lam*2)==lam*2:
                pass
                config_plot_xyz(filename, tag=r"$\lambda=%.1f,\bar{\kappa}=%.0f$" % (lam,karg),Format="png")
                #time.sleep(10)
                #config_plot_xyz(filename, r"$\tau_g^*=%.1f$" % tau0,fix_index=Fix_index)

        #print("sleeping...")
        #time.sleep(10)
    Os_pars_plot(foldername, pars, colors, alphas, "lam")
    Ieig_pars_plot(foldername, pars, colors, alphas, "lam")
    #Os_Nes_pars_plot(foldername, Ne1pars,Ne2pars,"lam")

    # additional test



if __name__ == '__main__':
    main()
