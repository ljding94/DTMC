import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
import scipy.fft
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from autocorrelation import *

def Os_pars_plot(foldername, pars, colors, alphas, mode):
    data,O_label, xLabel = Os_pars_data(foldername,"/O_",pars,mode)
    data = np.transpose(np.array(data), axes=(1, 0, 2))
    cpar, E_ave, E_tau, E_err = data[:4]
    Les_ave,Les_tau,Les_err = [],[],[]
    Ne = pars[0][1]
    for e in range(Ne):
        Les_ave.append(data[4+3*e])
        Les_tau.append(data[5+3*e])
        Les_err.append(data[6+3*e])
    Le_ave = np.sum(Les_ave,axis=0)
    Le_err = np.sqrt(np.sum(np.power(Les_err,2),axis=0))
    IdA_ave, IdA_tau, IdA_err,I2H_ave, I2H_tau, I2H_err,I2H2_ave, I2H2_tau, I2H2_err, Ikg_ave, Ikg_tau, Ikg_err, Ik2_ave, Ik2_tau, Ik2_err, Itau_ave, Itau_tau, Itau_err = data[7+3*(Ne-1):25+3*(Ne-1)]

    #Le_ave_N, I2H2_ave_N, Ik2_ave_N, Itau_ave_N, Itau2_ave_N = [], [], [], [], []
    #Le_err_N, I2H2_err_N, Ik2_err_N, Itau_err_N, Itau2_err_N = [], [], [], [], []

    F_ave,F_err=[],[]
    cpar_tsqN = []
    for i in range(len(pars)):
        if(mode=="L"):
            d0=1.2
            print(cpar[i])
            fa,fe = Chi2_gradient(d0*cpar[i],E_ave[i],E_err[i],4)
            F_ave.append(fa)
            F_err.append(fe)
    ppi = 72
    # LineWidth, FontSize, LabelSize = 1, 9, 8
    plt.figure()
    plt.rc('text', usetex=True)
    fig, axs = plt.subplots(10, 2, figsize=(
        246 / ppi*2, 246 / ppi * 5.5), sharex=True)  # , sharex=True
    #cpar_aj = cpar-np.outer([2.8, 2.0, 1.5, 0.8, 0], np.ones(len(cpar[0])))
    if(mode=="L"):
        d0=1.2
        O_cpar_plot(axs[-1,0], F_ave, F_err, O_label, "F", r"$\left<F\right>$",
                cpar, colors, alphas)
    O_cpar_plot(axs[0,0], E_ave, E_err, O_label, "E", r"$\left<E\right>$",
                cpar, colors, alphas)
    O_cpar_plot(axs[1,0], Le_ave, Le_err, O_label, "Le", r"$\left<\int ds\right>$",cpar, colors, alphas)
    if Ne==2:
        Le_ave_diff = np.abs(Les_ave[1]-Les_ave[0])
        Le_err_diff = np.sqrt(np.power(Les_err[1],2)+np.power(Les_err[0],2))
        O_cpar_plot(axs[1,1], Le_ave_diff, Le_err_diff, O_label, "Le_diff", r"$|\left<\int_0 ds\right>-\left<\int_1 ds\right>|$",cpar, colors, alphas)
        Ledif_ave,Ledif_tau,Ledif_err = data[25+3*(Ne-1):]
        O_cpar_plot(axs[2,1], Ledif_ave, Ledif_ave, O_label, "Le_diff'", r"$\left<|\int_0 ds-\int_1 ds|\right>$",cpar, colors, alphas)
        O_cpar_plot(axs[3,1], Le_ave_diff/Le_ave, Le_err_diff*0, O_label, "Le_diff/Le_ave", r"$|\left<\int_0 ds\right>-\left<\int_1 ds\right>|/sum$",cpar, colors, alphas)
        O_cpar_plot(axs[4,1], Ledif_ave/Le_ave, Ledif_ave*0, O_label, "Le_diff/Le_ave'", r"$\left<|\int_0 ds-\int_1 ds|\right>/sum$",cpar, colors, alphas)
        axs[4,1].set_xlabel(xLabel)

    O_cpar_plot(axs[2,0], IdA_ave, IdA_err, O_label, "IdA", r"$\left<\int dA\right>$",cpar, colors, alphas)
    O_cpar_plot(axs[3,0], I2H_ave, I2H_err, O_label, "I2H", r"$\left<\int dA (2H)\right>$",cpar, colors, alphas)

    O_cpar_plot(axs[4,0], I2H2_ave/(16*np.pi), I2H2_err/(16*np.pi), O_label, "I2H2", r"$\left<\int dA (2H)^2\right>/(16\pi)$",
                cpar, colors, alphas)
    O_cpar_plot(axs[5,0], Ikg_ave, Ikg_err, O_label, "Ikg", r"$\left<\int ds k_g\right>$",
                cpar, colors, alphas)
    O_cpar_plot(axs[6,0], Ik2_ave, Ik2_err, O_label, "Ik2", r"$\left<\int ds k^2\right>$",
                cpar, colors, alphas)
    O_cpar_plot(axs[7,0], Itau_ave, Itau_err, O_label, "Itau", r"$\left<\int ds \tau\right>$",
                cpar, colors, alphas)
    Itau_Le_err = Itau_ave / Le_ave * \
        np.sqrt(np.power(Itau_err / Itau_ave, 2) +
                np.power(Le_err / Le_ave, 2))
    O_cpar_plot(axs[8,0], Itau_ave / Le_ave, Itau_Le_err, O_label, "Itau_Le", r"$\left<\int ds \tau_g\right>/\left<\int ds\right>$",
                cpar, colors, alphas)
    #O_cpar_plot(axs[9,0], Itau2_ave, Itau2_err, O_label, "Itau2", r"$\left<\int ds \tau^2\right>$",cpar, colors, alphas)
    axs[8,0].set_xlabel(xLabel)
    #axs[0,0].xaxis.set_major_locator(MultipleLocator(2))
    axs[0,0].xaxis.set_minor_locator(MultipleLocator(1))
    lgd = axs[0,0].legend(loc="upper center",
                        bbox_to_anchor=(0.5, 0.25 * len(axs)))
    plt.tight_layout(pad=0)
    plt.savefig(foldername + "/Os_" + mode + ".pdf",
                format="pdf", bbox_extra_artists=(lgd,), bbox_inches='tight', transparent=True)
    plt.close()

    # individual plots, for slides
    '''
    ppi = 72
    # LineWidth, FontSize, LabelSize = 1, 9, 8
    plt.figure()
    plt.rc('text', usetex=True)
    fig, axs = plt.subplots(4, 2, figsize=(
        246 / ppi*2, 246 / ppi * 2))
    O_cpar_plot(axs[0,0], Le_ave, Le_err, O_label, "Le", r"$\left<\int ds\right>$",cpar, colors, alphas)
    if Ne==2:
        Le_ave_diff = np.abs(Les_ave[1]-Les_ave[0])
        Le_err_diff = np.sqrt(np.power(Les_err[1],2)+np.power(Les_err[0],2))
        O_cpar_plot(axs[0,1], Le_ave_diff, Le_err_diff, O_label, "Le_diff", r"$|\left<\int_0 ds\right>-\left<\int_1 ds\right>|$",cpar, colors, alphas)
        Ledif_ave,Ledif_tau,Ledif_err = data[25+3*(Ne-1):]
        O_cpar_plot(axs[1,1], Ledif_ave, Ledif_ave, O_label, "Le_diff'", r"$\left<|\int_0 ds-\int_1 ds|\right>$",cpar, colors, alphas)
        O_cpar_plot(axs[2,1], Le_ave_diff/Le_ave, Le_err_diff*0, O_label, "Le_diff/Le_ave", r"$|\left<\int_0 ds\right>-\left<\int_1 ds\right>|/sum$",cpar, colors, alphas)
        O_cpar_plot(axs[3,1], Ledif_ave/Le_ave, Ledif_ave*0, O_label, "Le_diff/Le_ave'", r"$\left<|\int_0 ds-\int_1 ds|\right>/sum$",cpar, colors, alphas)
        axs[3,1].set_xlabel(xLabel)

    O_cpar_plot(axs[1,0], IdA_ave, IdA_err, O_label, "IdA", r"$\left<\int dA\right>$",cpar, colors, alphas)
    O_cpar_plot(axs[2,0], I2H_ave, I2H_err, O_label, "I2H", r"$\left<\int dA (2H)\right>$",cpar, colors, alphas)

    O_cpar_plot(axs[3,0], I2H2_ave/(16*np.pi), I2H2_err/(16*np.pi), O_label, "I2H2", r"$\left<\int dA (2H)^2\right>/(16\pi)$",
                cpar, colors, alphas)
    axs[3,0].set_xlabel(xLabel)
    #axs[0,0].xaxis.set_major_locator(MultipleLocator(2))
    axs[0,0].xaxis.set_minor_locator(MultipleLocator(1))
    lgd = axs[0,0].legend(loc="upper center",
                        bbox_to_anchor=(0.5, 0.25 * len(axs)))
    plt.tight_layout(pad=0)
    plt.savefig(foldername + "/O_" + mode + ".pdf",
                format="pdf", bbox_extra_artists=(lgd,), bbox_inches='tight', transparent=True)
    plt.close()
    '''


def Ieig_pars_plot(foldername, pars, colors, alphas, mode):
    data,O_label, xLabel = Os_pars_data(foldername,"/Iij_",pars,mode)
    data = np.transpose(np.array(data), axes=(1, 0, 2))
    cpar, Ieig0_ave,Ieig0_tau,Ieig0_err,Ieig1_ave,Ieig1_tau,Ieig1_err,Ieig2_ave,Ieig2_tau,Ieig2_err = data
    print("cpar")
    print(cpar)
    print("Ieig0_ave")
    print(Ieig0_ave)
    ppi = 72
    # LineWidth, FontSize, LabelSize = 1, 9, 8
    plt.figure()
    plt.rc('text', usetex=True)
    fig, axs = plt.subplots(3, 1, figsize=(
        246 / ppi*1.5, 246 / ppi * 1.6),sharey=True,sharex=True)  # , sharex=True
    ylim_min = np.amin(Ieig0_ave-3*Ieig0_err)
    ylim_max = np.amax(Ieig2_ave+3*Ieig2_err)
    print("ylim_min,ylim_max")
    print(ylim_min,ylim_max)
    N  = np.array(pars)[:,0]
    for i in range(len(N)):
        Ieig0_ave[i]/=N[i]
        Ieig0_err[i]/=N[i]
        Ieig1_ave[i]/=N[i]
        Ieig1_err[i]/=N[i]
        Ieig2_ave[i]/=N[i]
        Ieig2_err[i]/=N[i]

    O_cpar_plot(axs[0], Ieig0_ave, Ieig0_err, O_label, "Ieig0", r"$\left<Ieig0\right>/N^2$",cpar, colors, alphas)
    O_cpar_plot(axs[1], Ieig1_ave, Ieig1_err, O_label, "Ieig1", r"$\left<Ieig1\right>/N^2$",cpar, colors, alphas)
    O_cpar_plot(axs[2], Ieig2_ave, Ieig2_err, O_label, "Ieig2", r"$\left<Ieig2\right>/N^2$",cpar, colors, alphas)
    axs[2].set_xlabel(xLabel)
    axs[0].set_ylim = (ylim_min,ylim_max)
    axs[0].xaxis.set_minor_locator(MultipleLocator(0.5))
    lgd = axs[0].legend(loc="upper center",bbox_to_anchor=(0.5, 1+0.15*len(axs)))
    #plt.tight_layout(pad=0)
    plt.savefig(foldername + "/Iij_" + mode + ".pdf")
    plt.close()


def Os_pars_data(foldername, head ,pars, mode):
    data, O_label = [], []
    if mode == "lam":
        xLabel = r"$\lambda\sigma_0/(k_B T)$"
        for i in range(len(pars)):
            N, Ne,L, kar,C0, karg,lam, B, Bc, tau0 = pars[i]
            filename = foldername + head+"MC_N%.0f_Ne%.0f_L%.0f_kar%.0f_C0%.1f_karg%.1f_lams_B%.1f_Bc%.0f_tau0%.2f_ana.txt" % (N, Ne, L, kar,C0, karg, B, Bc, tau0)
            print(filename)
            data.append(np.loadtxt(filename, skiprows=1,delimiter=",", unpack=True))
            O_label.append(r"$N=%.0f,Ne=%.0f,L=%.0f,\kappa=%.0f,C_0=%.1f,\bar{\kappa}=%.1f,B=%.0f,B'=%.0f,\tau_g^*=%.2f$" % (N, Ne,L, kar,C0, karg, B, Bc, tau0))
            #O_label.append(r"$\bar{\kappa}=%.1f$" % karg)
    elif mode == "tau0":
        xLabel = r"$\tau_g^*\sigma_0$"
        for i in range(len(pars)):
            N, Ne, L, kar,C0,karg,lam, B, Bc, tau0 = pars[i]
            filename = foldername + head+"MC_N%.0f_Ne%.0f_L%.0f_kar%.0f_C0%.1f_karg%.1f_lam%.1f_B%.1f_Bc%.0f_tau0s_ana.txt" % (N, Ne, L, kar, C0,lam, karg, B, Bc)
            data.append(np.loadtxt(filename, skiprows=1,delimiter=",", unpack=True))
            O_label.append(r"$N=%.0f, Ne=%.0f,L=%.0f,\kappa=%.0f,C_0=%.1f,\bar{\kappa}=%.1f,\lambda=%.1f,B=%.0f,B'=%.0f$" % (N, Ne, L, kar,C0, karg,lam, B, Bc))
            # E_parts_cpar_plot(data[i], pars[i], r"$\tau)g^*\sigma_0$", filename[:-4] + "_E.pdf")
    elif mode == "N":
        xLabel = r"$N$"
        for i in range(len(pars)):
            N, Ne, L, kar,C0, karg,lam, B, Bc, tau0 = pars[i]
            filename = foldername + head+"MC_Ns_Ne%.0f_L%.0f_kar%.0f_C0%.1f_karg%.1f_lam%.1f_B%.1f_Bc%.0f_tau0%.2f_ana.txt" % (Ne, L, kar, C0, karg,lam, B, Bc, tau0)
            data.append(np.loadtxt(filename, skiprows=1,delimiter=",", unpack=True))
            O_label.append(r"$L=%.0f,Ne=%.0f,\kappa=%.0f,C_0=%.1f,\bar{\kappa}=%.1f,\lambda=%.1f,B=%.0f,B'=%.0f,\tau_g^*=%.2f$" % (L,Ne, kar,C0, karg, lam, B, Bc, tau0))
    elif mode == "L":
        xLabel = r"$L$"
        for i in range(len(pars)):
            N, Ne, L, kar,C0,karg,lam, B, Bc, tau0 = pars[i]
            filename = foldername + head+"MC_N%.0f_Ne%.0f_Ls_kar%.0f_C0%.1f_karg%.1f_lam%.1f_B%.1f_Bc%.0f_tau0%.2f_ana.txt" % (N, Ne, kar, C0, karg,lam, B, Bc, tau0)
            data.append(np.loadtxt(filename, skiprows=1,delimiter=",", unpack=True))
            O_label.append(r"$N=%.0f,Ne=%.0f,\kappa=%.0f,C_0=%.1f,\bar{\kappa}=%.1f,\lambda=%.1f,B=%.0f,B'=%.0f,\tau_g^*=%.2f$" % (N, Ne, kar,C0, karg, lam, B, Bc, tau0))
    elif mode == "B":
        xLabel = r"$B/(K_BT\sigma_0)$"
        for i in range(len(pars)):
            N, Ne, L, kar, C0, karg,lam, B, Bc, tau0 = pars[i]
            filename = foldername + head+"MC_N%.0f_Ne%.0f_L%.0f_kar%.0f_C0%.1f_karg%.1f_lam%.1f_Bs_Bc%.0f_tau0%.2f_ana.txt" % (N, Ne, L, kar, C0, karg,lam, Bc, tau0)
            data.append(np.loadtxt(filename, skiprows=1,delimiter=",", unpack=True))
            O_label.append(r"$N=%.0f, Ne=%.0f,L=%.0f,\kappa=%.0f,C_0=%.1f,\bar{\kappa}=%.1f,\lambda=%.1f,B'=%.0f,\tau_g^*=%.2f$" % (N,Ne, L, kar, C0, karg, lam,Bc, tau0))
    elif mode == "Bc":
        xLabel = r"$B'/(K_BT\sigma_0)$"
        for i in range(len(pars)):
            N,Ne, L, kar,C0,  karg,lam, B, Bc, tau0 = pars[i]
            filename = foldername + head+"MC_N%.0f_Ne%.0f_L%.0f_kar%.0f_C0%.1f_karg%.1f_lam%.1f_B%.1f_Bcs_tau0%.2f_ana.txt" % (N, Ne, L, kar,C0, karg, lam, B, tau0)
            data.append(np.loadtxt(filename, skiprows=1,delimiter=",", unpack=True))
            O_label.append(r"$N=%.0f,Ne=%.0f,L=%.0f,\kappa=%.0f,C_0=%.1d,\bar{\kappa}=%.1f,\lambda=%.1f,B=%.0f,\tau_g^*=%.2f$" % (N, Ne, L, kar,C0, karg, lam, B, tau0))
    return data,O_label,xLabel


def O_cpar_plot(ax, O, O_err, O_label, O_name, O_latex, cpar, Colors, Alphas, Xscale="linear", Yscale="linear"):
    print("plotting %s vs cpar" % O_name)
    print(len(O),len(cpar),len(O_err))
    for i in range(len(O)):
        if(Colors):
            ax.errorbar(cpar[i], O[i], yerr=O_err[i], linestyle="None",
                        color=Colors[i], alpha=Alphas[i], label=O_label[i])
        else:
            ax.errorbar(cpar[i], O[i], yerr=O_err[i],
                        linestyle="-",marker=">", label=O_label[i])
    # ax.xlabel("$\lambda\sigma_0/k_B T$", fontsize=FontSize)
    ax.set_ylabel(O_latex)
    ax.set_xscale(Xscale)
    ax.set_yscale(Yscale)
    ax.tick_params(direction="in", top="on", right="on")


def Os_Nes_pars_plot(foldername, Ne1pars,Ne2pars, mode):
    Ne1data,Ne1O_label, Ne1xLabel = Os_pars_data(foldername,Ne1pars,mode)
    Ne2data,Ne2O_label, Ne2xLabel = Os_pars_data(foldername,Ne2pars,mode)
    Ne1data = np.transpose(np.array(Ne1data), axes=(1, 0, 2))
    Ne1cpar, Ne1E_ave, Ne1E_tau, Ne1E_err = Ne1data[:4]
    Ne1Les_ave,Ne1Les_tau,Ne1Les_err = [],[],[]
    for e in range(1):
        Ne1Les_ave.append(Ne1data[4+3*e])
        Ne1Les_tau.append(Ne1data[5+3*e])
        Ne1Les_err.append(Ne1data[6+3*e])
    Ne1Le_ave = np.sum(Ne1Les_ave,axis=0)
    Ne1Le_err = np.sqrt(np.sum(np.power(Ne1Les_err,2),axis=0))
    Ne2data = np.transpose(np.array(Ne2data), axes=(1, 0, 2))
    Ne2cpar, Ne2E_ave, Ne2E_tau, Ne2E_err = Ne2data[:4]
    Ne2Les_ave,Ne2Les_tau,Ne2Les_err = [],[],[]
    for e in range(2):
        Ne2Les_ave.append(Ne2data[4+3*e])
        Ne2Les_tau.append(Ne2data[5+3*e])
        Ne2Les_err.append(Ne2data[6+3*e])
    Ne2Le_ave = np.sum(Ne2Les_ave,axis=0)
    Ne2Le_err = np.sqrt(np.sum(np.power(Ne2Les_err,2),axis=0))

    fig_row = len(Ne1cpar)
    ppi = 72
    # LineWidth, FontSize, LabelSize = 1, 9, 8
    plt.figure()
    plt.rc('text', usetex=True)
    fig, axs = plt.subplots(fig_row, 2, figsize=(
        246 / ppi, 246 / ppi * fig_row*0.55))
    for i in range(len(axs)):
        E_ave = [Ne1E_ave[i],Ne2E_ave[i]]
        E_err = [Ne1E_err[i],Ne2E_err[i]]
        Le_ave = [Ne1Le_ave[i],Ne2Le_ave[i]]
        Le_err = [Ne1Le_err[i],Ne2Le_err[i]]
        O_label = [Ne1O_label[i],Ne2O_label[i]]
        cpar = [Ne1cpar[i],Ne2cpar[i]]
        O_cpar_plot(axs[i,0], E_ave, E_err, O_label, "E", r"$\left<E\right>$",cpar, None, None)
        O_cpar_plot(axs[i,1], Le_ave, Le_err, O_label, "Le", r"$\left<\int d s\right>$",cpar, None, None)
        axs[i,0].xaxis.set_minor_locator(MultipleLocator(1))
        lgd = axs[i,0].legend(loc="lower center")
    axs[-1,0].set_xlabel(Ne1xLabel)
    axs[-1,1].set_xlabel(Ne1xLabel)

    plt.tight_layout(pad=0)
    plt.savefig(foldername + "/Os_Nes_" + mode + ".pdf",
                format="pdf", bbox_extra_artists=(lgd,), bbox_inches='tight', transparent=True)
    plt.close()

def E_parts_cpar_plot(data, para, cpar_label, savename):
    N, kar, lam, karg, B, Bc, tau0 = para
    cpar, E_ave, E_tau, E_err, Le_ave, Le_tau, Le_err, I2H2_ave, I2H2_tau, I2H2_err, Ikg_ave, Ikg_tau, Ikg_err, Ik2_ave, Ik2_tau, Ik2_err, Itau_ave, Itau_tau, Itau_err, Itau2_ave, Itau2_tau, Itau2_err = data
    print("plotting energy by parts")
    ppi = 72
    plt.figure(figsize=(2 * 246 / ppi * 1, 246 / ppi * 0.8))
    plt.rc('text', usetex=True)
    plt.plot(cpar, -karg * Ikg_ave, "--",
             alpha=0.6, label=r"$-\bar{\kappa}\left<\int k_g ds\right>$")
    E_acum = -karg * Ikg_ave
    plt.fill_between(cpar, y1=E_acum, y2=lam * Le_ave + E_acum,
                     alpha=0.6, label=r"$\lambda \left<\int ds\right>$")
    E_acum += lam * Le_ave
    plt.fill_between(cpar, y1=E_acum, y2=0.5 * kar * I2H2_ave + E_acum,
                     alpha=0.6, label=r"$\frac{\kappa}{2} \left<\int (2H)^2 dA\right>$")
    E_acum += 0.5 * kar * I2H2_ave
    plt.fill_between(cpar, y1=E_acum, y2=0.5 * B * Ik2_ave + E_acum,
                     alpha=0.6, label=r"$\frac{B}{2}\left<\int  k^2 ds\right>$")
    E_acum += 0.5 * B * Ik2_ave
    plt.fill_between(cpar, y1=E_acum, y2=0.5 * Bc * (Itau2_ave - 2 * tau0 * Itau_ave + tau0 * tau0 * Le_ave) + E_acum,
                     alpha=0.6, label=r"$\frac{B'}{2}\left<\int (\tau-\tau_0)^2 ds\right>$")
    E_acum += 0.5 * Bc * (Itau2_ave - 2 * tau0 *
                          Itau_ave + tau0 * tau0 * Le_ave)
    plt.plot(cpar, E_ave, "--",
             alpha=0.6, label=r"$\left<E \right>$")
    plt.xlabel(cpar_label)
    plt.ylabel("$E_i$")
    plt.ylim(-50, 2 * N)
    plt.tick_params(direction="in", top="on", right="on")
    plt.title(r"$N=%.0f,\kappa=%.0f,\bar{\kappa}=%.1f,B=%.0f,B'=%.0f,\tau_g^*=%.2f$" % (
        N, kar, karg, B, Bc, tau0))
    plt.legend()
    plt.tight_layout(pad=0)
    plt.savefig(savename, format="pdf", transparant=True)
    plt.close()


def Cxx_stat_plot(foldername,Cmode, N_, Ne_, L_, kar_, lam_, karg_, B_, Bc_, tau0_, mode,p_int):

    Cxxs_es= []
    Cxxs_es_err= []
    Lens_es = []
    ave_xs_es = []
    ave_xs_es_err = []
    if (mode == "lam"):
        cpar = lam_
    elif (mode == "tau0"):
        cpar = tau0_
    elif (mode == "N"):
        cpar = N_
    elif (mode == "L"):
        cpar = L_
    elif (mode == "B"):
        cpar = B_
    elif (mode == "Bc"):
        cpar = Bc_
    for i in range(len(cpar)):
        if(mode == "lam"):
            N,Ne, L, kar, lam, karg, B, Bc, tau0 = N_, Ne_, L_, kar_, lam_[
                i], karg_, B_, Bc_, tau0_
        elif (mode == "tau0"):
            N,Ne, L, kar, lam, karg, B, Bc, tau0 = N_, Ne_, L_, kar_, lam_, karg_, B_, Bc_, tau0_[i]
        elif (mode == "N"):
            N, Ne,L, kar, lam, karg, B, Bc, tau0 = N_[
                i], Ne_,L_, kar_, lam_, karg_, B_, Bc_, tau0_
        elif (mode == "L"):
            N, Ne,L, kar, lam, karg, B, Bc, tau0 = N_, Ne_,L_[
                i], kar_, lam_, karg_, B_, Bc_, tau0_
        elif (mode == "B"):
            N,Ne, L, kar, lam, karg, B, Bc, tau0 = N_, Ne_,L_, kar_, lam_, karg_, B_[
                i], Bc_, tau0_
        elif (mode == "Bc"):
            N,Ne, L, kar, lam, karg, B, Bc, tau0 = N_, Ne_,L_, kar_, lam_, karg_, B_, Bc_[
                i], tau0_

        Cxxs_es.append([])
        Cxxs_es_err.append([])
        Lens_es.append([])
        ave_xs_es.append([])
        ave_xs_es_err.append([])

        for e in range(Ne) :
            print("autocorrelating ", N, Ne, L, kar, lam, karg, B, Bc, tau0)
            file2read = foldername + "/C"+Cmode+"_MC_N%.0f_Ne%.0f_L%.0f_kar%.0f_lam%.1f_karg%.1f_B%.1f_Bc%.0f_tau0%.2f_e%d.txt" % (
                N,Ne, L, kar, lam, karg, B, Bc, tau0,e)
            # len(edge_list)
            Len = np.genfromtxt(file2read, skip_header=1,
                                delimiter=",", usecols=0, dtype=int)
            Lens_es[i].append(np.average(Len))
            #
            # x = nc for ncnc
            # x = r for rr
            ave_x = np.genfromtxt(file2read, skip_header=1,
                                    delimiter=",", usecols=1, dtype=float)
            ave_xs_es[i].append(np.average(ave_x))
            rho, cov0 = autocorrelation_function_fft(ave_x)
            tau, tau_err = tau_int_cal_rho(rho,c=6)
            ave_xs_es_err[i].append(np.sqrt(2 * tau / len(ave_x) * cov0))

            # Cxx
            ucol={"rr":2,"ncnc":2}
            #Cxx = np.genfromtxt(file2read, skip_header=1,delimiter=",", usecols=range(min(Len)+2)[ucol[Cmode]:])
            Cxx = np.genfromtxt(file2read, skip_header=1,
                                delimiter=",", usecols=range(int(min(Len)/2)+1)[ucol[Cmode]:])
            Cxx_err = np.std(Cxx, axis=0)/np.sqrt(len(Cxx))
            Cxx = np.average(Cxx, axis=0)
            Cxxs_es[i].append(Cxx)
            Cxxs_es_err[i].append(Cxx_err)
    # for i in range(len(Cxxs)):
    # print(i, len(Cxxs[i]))
    if(mode == "lam"):
        savefile = foldername + "/C"+Cmode+"_MC_N%.0f_Ne%.0f_L%.0f_kar%.0f_lams_karg%.1f_B%.1f_Bc%.0f_tau0%.2f_plot.pdf" % (
            N,Ne, L, kar, karg, B, Bc, tau0)
    elif(mode == "tau0"):
        savefile = foldername + "/C"+Cmode+"_MC_N%.0f_Ne%.0f_L%.0f_kar%.0f_lam%.1f_karg%.1f_B%.1f_Bc%.0f_tau0s_plot.pdf" % (
            N,Ne, L, kar, lam, karg, B, Bc)
    elif(mode == "N"):
        savefile = foldername + \
            "/C"+Cmode+"_MC_Ns_L%.0f_kar%.0f_lam%.1f_karg%.1f_B%.1f_Bc%.0f_tau0%.2f_plot.pdf" % (
                Ne,L, kar, lam, karg, B, Bc, tau0)
    elif(mode == "L"):
        savefile = foldername + \
            "/C"+Cmode+"_MC_N%.0f_Ne%.0f_Ls_kar%.0f_lam%.1f_karg%.1f_B%.1f_Bc%.0f_tau0%.2f_plot.pdf" % (
                N,Ne, kar, lam, karg, B, Bc, tau0)
    elif(mode == "B"):
        savefile = foldername + \
            "/C"+Cmode+"_MC_N%.0f_Ne%.0f_L%.0f_kar%.0f_lam%.1f_karg%.1f_Bs_Bc%.0f_tau0%.2f_plot.pdf" % (
                N,Ne, L, kar, lam, karg, Bc, tau0)
    elif(mode == "Bc"):
        savefile = foldername + \
            "/C"+Cmode+"_MC_N%.0f_Ne%.0f_L%.0f_kar%.0f_lam%.1f_karg%.1f_B%.1f_Bcs_tau0%.2f_plot.pdf" % (
                N,Ne, L, kar, lam, karg, B, tau0)

    Cxxs_es = np.array(Cxxs_es)
    Cxxs_es_err = np.array(Cxxs_es_err)
    Lens_es = np.array(Lens_es)
    ave_xs_es = np.array(ave_xs_es)
    ave_xs_es_err = np.array(ave_xs_es_err)
    ppi = 72
    # LineWidth, FontSize, LabelSize = 1, 9, 8
    plt.figure()
    plt.rc('text', usetex=True)
    fig, axs = plt.subplots(max(Ne,2), 3, figsize=(
        246 / ppi*max(3,Ne), 246 / ppi * 1.5))  # , sharex=True
    print(np.shape(axs))
    for e in range(Ne):
        for i in range(len(Cxxs_es)):
            if(i % p_int == 0):
                if(Cmode=="ncnc"):
                    axs[e,0].errorbar(range(len(Cxxs_es[i,e]))/Lens_es[i,e], Cxxs_es[i,e],
                            yerr=Cxxs_es_err[i,e], label=mode+"=%.1f, edge%d" % (cpar[i],e))
                    axs[e,0].set_xlabel("s/L")
                elif(Cmode=="rr"):
                    Cuu = Cxxs_es[i,e]-np.power(ave_xs_es[i,e],2)
                    #axs[e,0].plot(range(len(Cuu))/Lens_es[i,e], np.ones(len(Cuu))*np.power(ave_xs_es[i],2),"--" , label=mode+"=%.1f, edge%d, <uu'>" % (cpar[i],e))
                    #print("Cxxs_es[i,e],i,cpar[i]",Cxxs_es[i,e],i,cpar[i])
                    axs[e,0].errorbar(range(len(Cxxs_es[i,e]))/Lens_es[i,e], Cuu,yerr=Cxxs_es_err[i,e], label=mode+"=%.1f, edge%d" % (cpar[i],e))
                    #axs[e,0].errorbar(range(len(Cxxs_es[i,e]))/ave_xs_es[i], Cuu,yerr=Cxxs_es_err[i,e], label=mode+"=%.1f, edge%d" % (cpar[i],e))
                    #axs[e,0].errorbar(range(len(Cxxs_es[i,e])), Cuu,yerr=Cxxs_es_err[i,e], label=mode+"=%.1f, edge%d" % (cpar[i],e))
                    axs[e,0].set_xlabel("s/L")
                    #axs[e,0].set_ylim(-0.7,1.6)
                    #axs[e,0].set_xlim(0.0,0.45)

                    uquq=scipy.fft.dct(Cuu)
                    #uquq=np.fft.rfft(Cuu)#[:int(len(Cuu)/2)-1]
                    #print("len(uquq),len(Cuu)",len(uquq),len(Cuu))
                    #freq=np.fft.rfftfreq(len(Cuu),d=1.0/Lens_es[i,e])#[:int(len(Cuu)/2)-1]
                    #freq=np.fft.rfftfreq(len(Cuu),d=1.0)#[:int(len(Cuu)/2)-1]
                    #freq=scipy.fft.rfftfreq(len(Cuu),d=1.0/Lens_es[i,e])#[:int(len(Cuu)/2)-1]
                    #print("len(Cuu),print(len(uquq)),print(len(freq))",len(Cuu),len(uquq),len(freq))
                    freq = np.array(range(len(uquq)))*np.pi/len(uquq)*Lens_es[i,e]
                    axs[e,1].plot(freq,uquq*np.power(freq,2),"+")
                    axs[e,2].plot(freq,uquq,"+")
                    axs[e,1].set_xlabel("qL")
                    axs[e,2].set_xlabel("qL")
                    axs[e,1].set_ylabel(r"$\left<u_q u_{-q}\right>(qL)^2$")
                    axs[e,2].set_ylabel(r"$\left<u_q u_{-q}\right>$")
                    axs[e,1].set_ylim(1,9e3)
                    axs[e,2].set_ylim(1e-4,30)
                    axs[e,1].set_xscale("log")
                    axs[e,1].set_yscale("log")
                    axs[e,2].set_xscale("log")
                    axs[e,2].set_yscale("log")
            # do the fourier transform

        #axs[e,0].set_ylim(-0.1,1)
        labelx = {"rr":r"$\left<\left<u(i)\cdot u(i+s)\right>_i\right>$","ncnc":r"$\left<\left<n_c(i)\cdot n_c(i+s)\right>_i\right>$"}
        axs[e,0].set_ylabel(labelx[Cmode])
        # plt.yscale("log")
        axs[e,0].legend()
        print("np.shape(cpar)",np.shape(cpar))
        print("np.shape(ave_xs_es[:,e])",np.shape(ave_xs_es[:,e]))
        if(Cmode=="ncnc"):
            print(ave_xs_es[:,e])
            axs[e,1].errorbar(cpar,ave_xs_es[:,e],yerr=ave_xs_es_err[:,e])
            axs[e,1].set_xlabel(mode)
            axs[e,1].set_ylim(0,1)
            axs[e,1].set_ylabel(r"$\left<|\frac{\sum_{\partial\mathcal{M} \hat{n}_c }}{|\partial\mathcal{M}} |\right>$")
            axs[e,2].remove()
        if Ne==1:
            for x in range(len(axs[1])):
                axs[1,x].remove()

    plt.savefig(savefile, format="pdf", bbox_inches='tight', transparent=True)
    plt.close()


def Chi2_gradient(x,y,yerr,k):

    #input x,y (n,1) array
    #output b,berr (n-2,1) array
    #fit near 2k+1 point using chi2 method

    b = np.zeros(len(x))
    berr = np.zeros(len(x))

    # first k points

    for i in range(k):
        x_t = np.array(x[:i+k+1])
        y_t = np.array(y[:i+k+1])
        sig_t = np.array(yerr[:i+k+1])
        Delta = np.sum(1/sig_t/sig_t)*np.sum(np.power(x_t/sig_t,2))-np.power(np.sum(x_t/sig_t/sig_t),2)
        b[i]=1/Delta*(np.sum(np.power(sig_t,-2))*np.sum(x_t*y_t/sig_t/sig_t)-np.sum(x_t/sig_t/sig_t)*np.sum(y_t/sig_t/sig_t))
        berr[i]=np.sqrt(1/Delta*np.sum(np.power(sig_t,-2)))


    for i in range(len(x)-k)[k:]:
        x_t = np.array(x[i-k:i+k+1])
        y_t = np.array(y[i-k:i+k+1])
        sig_t = np.array(yerr[i-k:i+k+1])
        Delta = np.sum(np.power(sig_t,-2))*np.sum(np.power(x_t/sig_t,2))-np.power(np.sum(x_t/sig_t/sig_t),2)
        #print("Delta: ",Delta)
        b[i]=1/Delta*(np.sum(1/sig_t/sig_t)*np.sum(x_t*y_t/sig_t/sig_t)-np.sum(x_t/sig_t/sig_t)*np.sum(y_t/sig_t/sig_t))
        #print("b[%.0f]="%i,b[i])
        berr[i]=np.sqrt(1/Delta*np.sum(np.power(sig_t,-2)))


    # last k points
    for i in range(k):
        i=-i-1
        x_t = np.array(x[i-k:])
        y_t = np.array(y[i-k:])
        sig_t = np.array(yerr[i-k:])
        Delta = np.sum(1/sig_t/sig_t)*np.sum(np.power(x_t/sig_t,2))-np.power(np.sum(x_t/sig_t/sig_t),2)
        b[i]=1/Delta*(np.sum(np.power(sig_t,-2))*np.sum(x_t*y_t/sig_t/sig_t)-np.sum(x_t/sig_t/sig_t)*np.sum(y_t/sig_t/sig_t))
        #print("b[%.0f]="%i,b[i])
        berr[i]=np.sqrt(1/Delta*np.sum(np.power(sig_t,-2)))

    return b,berr

