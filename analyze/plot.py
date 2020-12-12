import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib.patches import Circle


def config_plot_xyz(filename,tag="", Format="pdf",lim=15,fix_index=None):
    print("plotting",filename)
    data = np.loadtxt(filename, skiprows=19, delimiter=",", unpack=True)
    x, y, z, enum, en0, en1 = data[7:13]
    x_min, x_max = np.min(x),np.max(x)
    y_min, y_max = np.min(y),np.max(y)
    z_min, z_max = np.min(z),np.max(z)
    alpha_xy = 0.9*(z-z_min+0.1)/(z_max-z_min+0.1)+0.1
    alpha_zx = 0.9*(y-y_min+0.1)/(y_max-y_min+0.1)+0.1
    ns = np.transpose(data[12:])
    ens = np.array([en0, en1])
    fig = plt.figure(figsize=(10, 5))
    #ax_xy = fig.add_subplot(111, aspect="equal")
    ax_xy = fig.add_subplot(121, aspect="equal")
    ax_zx = fig.add_subplot(122, aspect="equal")
    for i in range(len(ns)):
        for j in range(len(ns[0])):
            if ns[i, j] != -1:
                ax_xy.plot([x[i],x[int(ns[i, j])]], [y[i],y[int(ns[i, j])]], color="tomato",alpha=alpha_xy[i])
                ax_zx.plot([z[i],z[int(ns[i, j])]], [x[i],x[int(ns[i, j])]], color="tomato",alpha=alpha_xy[i])

            #elif ns[i, j] == -1 and j<min_j:
            #    min_j=j
            #    print("bead#",i,ns[i])
    ecolors = ["blue","purple","green"]
    for i in range(len(ens)):
        for j in range(len(en0)):
            if ens[i, j] != -1:
                ax_xy.plot([x[j], x[int(ens[i, j])]], [
                    y[j], y[int(ens[i, j])]], "-", color=ecolors[int(enum[j])], alpha=alpha_xy[j])
                ax_zx.plot([z[j], z[int(ens[i, j])]],[x[j], x[int(ens[i, j])]], "-", color=ecolors[int(enum[j])], alpha=alpha_zx[j])
    #plot fixed bead (if have)
    if fix_index:
        ax_xy.plot([x[fix_index[0]], x[fix_index[1]]], [y[fix_index[0]], y[fix_index[1]]], marker="o",linestyle="None", color="purple")
        ax_zx.plot([z[fix_index[0]], z[fix_index[1]]], [x[fix_index[0]], x[fix_index[1]]], marker="o",linestyle="None", color="purple")
    '''
    for i in range(len(x)):
        ax_xy.annotate(i, (x[i], y[i]),fontsize=5)
        ax_zx.annotate(i, (z[i], x[i]),fontsize=5)
    '''
    '''
    r = 0.5 * np.ones(len(x))
    for j in range(len(x)):
        xc, yc, zc, rc = x[j], y[j], z[j], r[j]
        ax_xy.add_artist(Circle(xy=(xc, yc), radius=rc, alpha=alpha_xy[j]))
        ax_zx.add_artist(Circle(xy=(zc, xc), radius=rc, alpha=alpha_zx[j]))
    '''
    ax_xy.set_xlim(x_min-2, x_max+2)
    ax_xy.set_ylim(y_min-2, y_max+2)
    #ax_xy.set_aspect("equal")
    #ax_xy.set_yticks([])
    #ax_xy.set_xticks([])
    #ax_xy.set_frame_on(False)
    ax_xy.set_title("XY  "+tag, fontsize=25)
    ax_zx.set_xlim(z_min-2, z_max+2)
    ax_zx.set_ylim(x_min-2, x_max+2)
    ax_zx.set_title("ZX")

    # plt.savefig(filename[:-4] + "_xy.pdf", dpi=300, format="pdf")
    plt.savefig(filename[:-4] + "_xyz."+Format, dpi=50,
                format=Format, bbox_inches='tight',transparent=False)
    plt.close()

def config_plot_xyz_seq(filename,Seq):
    for i in range(Seq):
        config_plot_xyz(filename[:-4]+"_%d.txt"%i,Format="png")


def config_plot3D(filename):
    data = np.loadtxt(filename, skiprows=19, delimiter=",", unpack=True)
    x, y, z, enum,en0, en1 = data[6:12]
    ns = data[12:]
    ens = np.array([en0, en1])
    fig = plt.figure(figsize=(5, 5))
    ax = plt.axes(projection="3d")
    for i in range(len(ns)):
        for j in range(len(ns[0])):
            if ns[i, j] != -1:
                ax.plot3D([x[j], x[int(ns[i, j])]], [y[j], y[int(ns[i, j])]], [z[j], z[int(ns[i, j])]], "-",  color="tomato")
                #ax.quiver(x[j], y[j], z[j], 0.5 * (x[int(ns[i, j])] - x[j]), 0.5 * (y[int(ns[i, j])] - y[j]), 0.5 * (z[int(ns[i, j])] - z[j]), color="tomato")
    ecolors = ["blue","purple","green"]
    for i in range(len(ens)):
        for j in range(len(en0)):
            if ens[i, j] != -1:
                ax.plot3D([x[j], x[int(ens[i, j])]], [
                    y[j], y[int(ens[i, j])]], [
                    z[j], z[int(ens[i, j])]], "-", color=ecolors[int(enum[j])], alpha=0.7)
    #ax.set_xlim(-t*max_xy, t*max_xy)
    #ax.set_ylim(-t*max_xy, t*max_xy)
    #ax.set_zlim(-t*max_xy, t*max_xy)
    ax.scatter3D(x, y, z, s=[1 for i in range(len(x))])
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    plt.show()
    # plt.savefig(filename[:-4] + "_3D.png", dpi=300)
    plt.close()

def plot_sketch():
    x1 = np.linspace(0, 1, 100)
    y1 = np.sqrt(1 - np.power(x1, 2))
    x2 = np.linspace(1 / np.sqrt(2), 2, 100)
    y2 = np.linspace(1 / np.sqrt(2), 2, 100)
    plt.figure(figsize=(10, 5))
    # plt.plot(x1, y1, "-", color="tomato")
    # plt.plot(x2, y2, "-", color="royalblue")
    # plt.tick_params(axis='both', left='off', top='off',right = 'off', bottom = 'off')
    plt.xticks((4, 8))
    plt.yticks((4, 6, 8))
    plt.xlim(1.9, 10.1)
    plt.ylim(1.9, 10.1)
    plt.xlabel(r"$\kappa/(k_B T) $", fontsize=18)
    plt.ylabel(r"$\lambda\sigma_0/(k_B T)$", fontsize=18)
    plt.tick_params(direction="out", top="on", right="on", labelsize=15)
    # plt.annotate("jagged", (0.3, 0.3), fontsize=18)
    # plt.annotate("sphere", (0.3, 1.3), fontsize=18)
    # plt.annotate("disk", (1.3, 0.3), fontsize=18)
    plt.savefig("sketch_phase.pdf", dpi=300, format="pdf", bbox_inches='tight')
    # plt.show()
    plt.close()

def O_kappas_lams_color_plot(data, savefilename):
    kappa, lam, L_ave, L_std, En_ave, En_std = data
    kappa_sort = sorted(set(kappa))
    print(kappa_sort)
    lam_sort = sorted(set(kappa))
    kappa_dis = kappa_sort[1] - kappa_sort[0]
    kappa_min, kappa_max = kappa_sort[0], kappa_sort[-1]
    lam_dis = lam_sort[1] - lam_sort[0]
    lam_min, lam_max = lam_sort[0], lam_sort[-1]
    L_ave_grid = np.empty((len(kappa_sort), len(lam_sort)))
    L_ave_grid.fill(np.nan)

    for i in range(len(kappa)):
        kappa_pos = int((kappa[i] - kappa_min) / kappa_dis)
        lam_pos = int((lam[i] - lam_min) / lam_dis)
        L_ave_grid[lam_pos][kappa_pos] = L_ave[i]
    plt.figure()
    kappa_plot = np.append(kappa_sort - 0.5 * kappa_dis,
                           kappa_sort[-1] + 0.5 * kappa_dis)
    lam_plot = np.append(lam_sort - 0.5 * lam_dis,
                         lam_sort[-1] + 0.5 * lam_dis)
    plt.pcolormesh(kappa_plot, lam_plot, L_ave_grid)
    # plt.xticks(kappa_sort)
    # plt.yticks(lam_sort)
    plt.xlabel(r"$\kappa/(k_B T)$", fontsize=18)
    plt.ylabel(r"$\lambda\sigma_0/(k_B T)$", fontsize=18)
    plt.tick_params(direction="in", top="on", right="on", labelsize=15)
    cb = plt.colorbar()
    cb.ax.set_title(r"$L/\sigma_0$", fontsize=18)
    plt.show()
    #plt.savefig(savefilename, dpi=300, format="pdf", bbox_inches='tight')
    plt.close()

def O_lams_line_multi_kappa_plot(filenames, kappas, savefilename):

    fig = plt.figure(figsize=(5, 10))
    ax = fig.add_subplot(311)
    ax1 = fig.add_subplot(312, sharex=ax)
    ax2 = fig.add_subplot(313, sharex=ax)
    for i in range(len(filenames)):
        data = np.loadtxt(filenames[i], skiprows=1, delimiter=",", unpack=True)
        lams, L_aves, L_stds, i2H_aves, i2H_stds, En_aves, En_stds = data
        # lams,L_ave,L_std,i2H_ave,i2H_std,En_ave,En_std
        ax.errorbar(lams - kappas[i], L_aves, yerr=L_stds,
                    label=r"$\kappa/(k_B T)=%g$" % kappas[i])
        ax1.errorbar(lams - kappas[i], i2H_aves, yerr=i2H_stds,
                     label=r"$\kappa/(k_B T)=%g$" % kappas[i])
        ax2.errorbar(lams - kappas[i], En_aves, yerr=En_stds,
                     label=r"$\kappa/(k_B T)=%g$" % kappas[i])
    #ax.set_xlabel(r"$(\lambda\sigma_0-\kappa)/(k_B T)$")
    ax.set_ylabel(r"$L/\sigma_0$")
    ax.legend()
    ax1.set_ylabel(r"$\int (2H)^2 dA$")
    ax1.legend()
    ax2.set_ylabel(r"$E$")
    ax2.set_xlabel(r"$(\lambda\sigma_0-\kappa)/(k_B T)$")
    ax2.legend()
    # plt.show()
    plt.savefig(savefilename, dpi=300)
    plt.close()


def O_kar_lam_MCstep_plot(Ls, i2Hs, Ens, savefile):
    MCstep = np.arange(len(Ls))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(MCstep, Ls, "d", label=r"$L\sigma_0$")
    ax.plot(MCstep, i2Hs, "s", label=r"$\int (2H)^2 dA$")
    ax.plot(MCstep, Ens, "s", label=r"$E/k_BT$")
    ax.set_ylim(0,500)
    ax.set_xlabel("MC steps")
    plt.legend()
    plt.savefig(savefile,format="pdf",transparent=True)
    plt.close()

def O_kappa_En_dis_lam_plot(Ens,lams,N,savefile):
    dlamN = int(len(lams)/N)
    plt.figure(figsize=(4,2*dlamN))
    for i in range(dlamN):
        plt.subplot(dlamN,1,i+1)
        Ens_p = Ens[i::dlamN]
        lams_p = lams[i::dlamN]
        for j in range(len(Ens_p)):
            plt.hist(Ens_p[j],bins=100,histtype="step",label=r"$\lambda=%.1f$"%lams_p[j])
        plt.legend()
    plt.xlabel(r"$E/k_B T$")
    plt.ylabel("distribution")
    plt.tight_layout()
    plt.savefig(savefile)
    plt.close()



def O_scatter_kappa_lam_plot(filenames, kappas, N, savefile, filenames_rev=None, kappas_rev=None):
    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(111)
    for i in range(len(filenames)):
        data = np.loadtxt(filenames[i], skiprows=1, delimiter=",", unpack=True)
        lams, L_aves, L_stds, i2H_aves, i2H_stds, En_aves, En_stds = data
        plt.scatter(kappas[i] * np.ones(len(lams)),lams,
                    s=10 * np.ones(len(lams)), c=L_aves, marker="s", alpha=0.6)
        plt.clim(0, N)
    if(filenames_rev):
        for i in range(len(filenames_rev)):
            data = np.loadtxt(filenames_rev[i], skiprows=1,
                              delimiter=",", unpack=True)
            lams, L_aves, L_stds, i2H_aves, i2H_stds, En_aves, En_stds = data
            plt.scatter(kappas_rev[i] * np.ones(len(lams)) + 0.2,
                        lams, s=10 * np.ones(len(lams)), c=L_aves, marker="s", alpha=0.6)
            plt.clim(0, N)

    xp = np.linspace(5, 14, 100)
    plt.plot(xp, 4 * np.sqrt(np.pi/N) * xp, "--",
             label=r"$\frac{4 \sqrt{\pi}}{\sqrt{N}}\kappa/(k_B T)$")
    plt.xlabel(r"$\kappa/(k_B T)$")
    plt.xlim(0, 14)
    plt.ylabel(r"$\lambda\sigma_0/(k_B T)$")
    plt.ylim(0, 16)
    cbar = plt.colorbar()
    cbar.set_ticks([0,0.2*N,0.4*N,0.6*N,0.8*N,N])
    cbar.set_ticklabels(["0",r"$0.2N\sigma_0$",r"$0.4N\sigma_0$",r"$0.6N\sigma_0$",r"$0.8N\sigma_0$",r"$N\sigma_0$"])
    plt.legend()
    plt.tight_layout()
    #plt.show()
    plt.savefig(savefile, format="pdf",transparent=True)
    plt.close()


def O_scatter_kappa_lam_ave_plot(folders, filenames, kappas, N, savefile):
    fig = plt.figure(figsize=(4, 4))
    ax = fig.add_subplot(111)
    #ax1 = fig.add_subplot(312, sharex=ax)
    #ax2 = fig.add_subplot(313, sharex=ax)
    for i in range(len(filenames)):
        L_aves, i2H_aves, En_aves = [], [], []
        for j in range(len(folders)):
            filename = folders[j] + "/" + filenames[i]
            data = np.loadtxt(filename, skiprows=1, delimiter=",", unpack=True)
            lams, L_aves_, L_stds, i2H_aves_, i2H_stds, En_aves_, En_stds_ = data
            L_aves.append(L_aves_)
            i2H_aves.append(i2H_aves_)
            En_aves.append(En_aves_)
        L_aves = np.average(L_aves, axis=0)
        i2H_aves = np.average(i2H_aves, axis=0)
        En_aves = np.average(En_aves, axis=0)
        plt.scatter(kappas[i] * np.ones(len(lams)), np.sqrt(N) * lams,
                    c=L_aves, label="N=%d" % N)
        plt.clim(0, N)
        xp = np.linspace(5, 13, 100)
        plt.plot(xp, 4 * np.sqrt(np.pi) * xp, "--")
        #ax1.errorbar(lams - kappas[i], i2H_aves, yerr=i2H_stds,label=r"$\kappa/(k_B T)=%g$" % kappas[i])
        #ax2.errorbar(lams - kappas[i], En_aves, yerr=En_stds,label=r"$\kappa/(k_B T)=%g$" % kappas[i])
    plt.xlabel(r"$\kappa/(k_B T)$")
    plt.ylabel(r"$\sqrt{N}\lambda\sigma_0/(k_B T)$")
    plt.colorbar()
    plt.legend()
    #ax1.set_ylabel(r"$\int (2H)^2 dA$")
    # ax1.legend()
    # ax2.set_ylabel(r"$E$")
    #ax2.set_xlabel(r"$(\lambda\sigma_0-\kappa)/(k_B T)$")
    # ax2.legend()
    # plt.show()
    plt.savefig(savefile, dpi=300)
    plt.close()

def O_L_plot(filename,d0):
    L,En_ave,En_std,L_ave,L_std=np.loadtxt(filename,skiprows=1,delimiter=",",unpack=True)
    Fn_ave = [(En_ave[i+1]-En_ave[i])/d0 for i in range(len(En_ave)-1)]
    print(Fn_ave)
    Fn_ave=np.gradient(En_ave,1)
    print(Fn_ave)
    #Ln_F = L[:-1]+0.5
    fig = plt.figure(figsize=(5,8))
    ax_E = fig.add_subplot(211)
    ax_F = fig.add_subplot(212)
    ax_E.errorbar(L*d0,En_ave,yerr=En_std,label=r"$E/k_B T$")
    ax_E.errorbar(L*d0,L_ave,yerr=L_std,label=r"$L_e/\sigma_0$")
    ax_E.set_xlabel(r"$l_{fix}/\sigma_0$")
    ax_F.plot(L*d0,Fn_ave,label=r"$F\sigma_0/k_B T$")
    ax_F.set_xlabel(r"$l_{fix}/\sigma_0$")
    plt.legend(frameon=0)
    plt.savefig(filename[:-4]+".pdf",format="pdf",transparent=True)
    plt.show()

def Shape_helix_energy():
    y0 = 5*1.2*np.sqrt(3)/2
    q = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]
    q = np.array(q)
    E_m = [0.359522,1.24792,2.33724,3.40905,4.33224,5.00607,5.35523,5.34456]
    E_m_nci = [0.352461,1.25733,2.4577,3.80575,5.2411,6.7206,8.19097,9.5866]
    #E_m2 = [1.43809,4.99168,9.34896,13.6362,17.3289]
    E_m = np.array(E_m)
    E_m_nci = np.array(E_m_nci)
    #E_m2 = np.array(E_m2)/4
    qs = np.linspace(0,0.8,100)
    #E_p = 25*np.power(qs,2)/np.sqrt(27/4*np.power(qs,2)+1)
    E_p = 25*qs/(2*np.sqrt(2)*3*1.2)
    plt.figure()
    #plt.plot(q,E_m,"d",label="measured1")
    #plt.plot(q,E_m_nci,"d",label="measured nci")
    plt.plot(q,np.sqrt(E_m_nci/(25*1.2)),"-",label="theory")
    plt.plot(q,q,"-",label="theory")
    plt.xlabel("q*d0")
    plt.ylabel("integral")
    plt.legend()
    plt.show()
    #plt.savefig("torsion_test.pdf",format="pdf",transparent=True)
    plt.close()

def autocorrelation_plot(rho,tau_int,savefile):
    t = np.linspace(0,1000,1000)
    plt.figure()
    plt.plot(range(1000),rho[:1000],"d")
    plt.plot(t,np.exp(-t/tau_int),"--")
    plt.savefig(savefile,format="pdf",transparent=True)
    plt.close()
def O_scatter_lam_tau0_plot(filenames,N,L, kar, lam,savefile):
    fig = plt.figure(figsize=(12, 8))
    #ax = fig.add_subplot(111)
    axE = fig.add_subplot(231)
    axLe = fig.add_subplot(232)
    axI2H2 = fig.add_subplot(233)
    axIkg = fig.add_subplot(234)
    axIk2 = fig.add_subplot(235)
    axItau = fig.add_subplot(236)
    for i in range(len(filenames)):
        data = np.loadtxt(filenames[i], skiprows=1, delimiter=",", unpack=True)
        tau0,E_ave,E_tau,E_err,Le_ave,Le_tau,Le_err,I2H2_ave,I2H2_tau,I2H2_err,Ikg_ave,Ikg_tau,Ikg_err,Ik2_ave,Ik2_tau,Ik2_err,Itau_ave,Itau_tau,Itau_err = data
        Esc=axE.scatter(lam[i] * np.ones(len(tau0)),tau0,
                    s=10 * np.ones(len(tau0)), c=E_ave, marker="s", alpha=0.6)
        Lesc=axLe.scatter(lam[i] * np.ones(len(tau0)),tau0,
                    s=10 * np.ones(len(tau0)), c=Le_ave,marker="s", alpha=0.6)
        I2H2sc=axI2H2.scatter(lam[i] * np.ones(len(tau0)),tau0,
                    s=10 * np.ones(len(tau0)), c=I2H2_ave,vmin=-np.pi,vmax=5*np.pi, marker="s", alpha=0.6)
        #Ikgsc=axIkg.scatter(lam[i] * np.ones(len(tau0)),tau0,
        #            s=10 * np.ones(len(tau0)), c=Ikg_ave,vmin=-2*np.pi,vmax=4*np.pi, marker="s", alpha=0.6)
        Ikgsc=axIkg.scatter(lam[i] * np.ones(len(tau0)),tau0,
                    s=10 * np.ones(len(tau0)), c=Ikg_ave, marker="s", alpha=0.6)
        Ik2sc=axIk2.scatter(lam[i] * np.ones(len(tau0)),tau0,
                    s=10 * np.ones(len(tau0)), c=Ik2_ave, marker="s", alpha=0.6)
        #Itausc=axItau.scatter(lam[i] * np.ones(len(tau0)),tau0,
        #            s=10 * np.ones(len(tau0)), c=Itau_ave-tau0*Le_ave , marker="s", alpha=0.6)
        Itausc=axItau.scatter(lam[i] * np.ones(len(tau0)),tau0,
                    s=10 * np.ones(len(tau0)), c=Itau_ave/Le_ave , marker="s", alpha=0.6)
    axE.set_title(r"$E$")
    axLe.set_title(r"$\int d s$")
    axI2H2.set_title(r"$\int d A (2H)^2$")
    axIkg.set_title(r"$\int d s k_g$")
    axIk2.set_title(r"$\int d s k^2$")
    #axItau.set_title(r"$\int d s (\tau-\tau_0)$")
    axItau.set_title(r"$\int d s \tau/\int d s$")
    Ecbar=fig.colorbar(Esc, ax=axE)
    Lecbar=fig.colorbar(Lesc, ax=axLe)
    Lecbar.set_ticks([0,0.2*N,0.4*N,0.6*N,0.8*N,N])
    #Lecbar.set_ticklabels(["0",r"$0.2N\sigma_0$",r"$0.4N\sigma_0$",r"$0.6N\sigma_0$",r"$0.8N\sigma_0$",r"$N\sigma_0$"])
    I2H2cbar=fig.colorbar(I2H2sc, ax=axI2H2)
    I2H2cbar.set_ticks([0,np.pi,2*np.pi,3*np.pi,4*np.pi])
    I2H2cbar.set_ticklabels([r"$-\pi$","0",r"$\pi$",r"$2\pi$",r"$3\pi$",r"$4\pi$",r"$5\pi$"])
    Ikgcbar=fig.colorbar(Ikgsc, ax=axIkg)
    #Ikgcbar.set_ticks([-2*np.pi,-np.pi,0,np.pi,2*np.pi,3*np.pi,4*np.pi])
    #Ikgcbar.set_ticklabels([r"$-2\pi$",r"$-\pi$","0",r"$\pi$",r"$2\pi$",r"$3\pi$",r"$4\pi$"])
    Ik2cbar=fig.colorbar(Ik2sc, ax=axIk2)
    Itaucbar=fig.colorbar(Itausc, ax=axItau)

    #ax.set_xlabel(r"$\lambda\sigma_0/(k_B T)$")
    #ax.set_ylabel(r"$\tau_0\sigma_0$")
    fig.text(0.5,0.04,r"$\lambda\sigma_0/(k_B T)$", ha="center")
    fig.text(0.04,0.5,r"$\tau_0\sigma_0$",va="center",rotation="vertical")
    #plt.tight_layout()
    #plt.show()
    plt.savefig(savefile, format="pdf",transparent=True)
    plt.close()

def O_scatter_Bs_Bcs_plot(filenames,N,Bs,savefile):
    fig = plt.figure(figsize=(12, 8))
    #ax = fig.add_subplot(111)
    axE = fig.add_subplot(231)
    axLe = fig.add_subplot(232)
    axI2H2 = fig.add_subplot(233)
    axIkg = fig.add_subplot(234)
    axIk2 = fig.add_subplot(235)
    axItau = fig.add_subplot(236)
    for i in range(len(filenames)):
        data = np.loadtxt(filenames[i], skiprows=1, delimiter=",", unpack=True)
        Bcs, E_ave, E_tau, E_err, Le_ave, Le_tau, Le_err, I2H2_ave, I2H2_tau, I2H2_err, Ikg_ave, Ikg_tau, Ikg_err, Ik2_ave, Ik2_tau, Ik2_err, Itau_ave, Itau_tau, Itau_err, Itau2_ave, Itau2_tau, Itau2_err, Tau_ave, Tau_tau, Tau_err = data
        print("Bs",Bs[i] * np.ones(len(Bcs)))
        print("Bcs",Bcs)
        print(E_ave)
        Esc=axE.scatter(Bs[i] * np.ones(len(Bcs)),Bcs,
                    s=10 * np.ones(len(Bcs)), c=E_ave, marker="s", alpha=0.6)
        Lesc=axLe.scatter(Bs[i] * np.ones(len(Bcs)),Bcs,
                    s=10 * np.ones(len(Bcs)), c=Le_ave,marker="s", alpha=0.6)
        I2H2sc=axI2H2.scatter(Bs[i] * np.ones(len(Bcs)),Bcs,
                    s=10 * np.ones(len(Bcs)), c=I2H2_ave,vmin=-np.pi,vmax=5*np.pi, marker="s", alpha=0.6)
        #Ikgsc=axIkg.scatter(Bs[i] * np.ones(len(Bcs)),Bcs,
        #            s=10 * np.ones(len(Bcs)), c=Ikg_ave,vmin=-2*np.pi,vmax=4*np.pi, marker="s", alpha=0.6)
        Ikgsc=axIkg.scatter(Bs[i] * np.ones(len(Bcs)),Bcs,
                    s=10 * np.ones(len(Bcs)), c=Ikg_ave, marker="s", alpha=0.6)
        Ik2sc=axIk2.scatter(Bs[i] * np.ones(len(Bcs)),Bcs,
                    s=10 * np.ones(len(Bcs)), c=Ik2_ave, marker="s", alpha=0.6)
        #Itausc=axItau.scatter(Bs[i] * np.ones(len(Bcs)),Bcs,
        #            s=10 * np.ones(len(Bcs)), c=Itau_ave-Bcs*Le_ave , marker="s", alpha=0.6)
        Itausc=axItau.scatter(Bs[i] * np.ones(len(Bcs)),Bcs,
                    s=10 * np.ones(len(Bcs)), c=Itau_ave/Le_ave , marker="s", alpha=0.6)
    axE.set_title(r"$E$")
    axLe.set_title(r"$\int d s$")
    axI2H2.set_title(r"$\int d A (2H)^2$")
    axIkg.set_title(r"$\int d s k_g$")
    axIk2.set_title(r"$\int d s k^2$")
    #axItau.set_title(r"$\int d s (\tau-\tau_0)$")
    axItau.set_title(r"$\int d s \tau/\int d s$")
    Ecbar=fig.colorbar(Esc, ax=axE)
    Lecbar=fig.colorbar(Lesc, ax=axLe)
    Lecbar.set_ticks([0,0.2*N,0.4*N,0.6*N,0.8*N,N])
    #Lecbar.set_ticklabels(["0",r"$0.2N\sigma_0$",r"$0.4N\sigma_0$",r"$0.6N\sigma_0$",r"$0.8N\sigma_0$",r"$N\sigma_0$"])
    I2H2cbar=fig.colorbar(I2H2sc, ax=axI2H2)
    I2H2cbar.set_ticks([0,np.pi,2*np.pi,3*np.pi,4*np.pi])
    I2H2cbar.set_ticklabels([r"$-\pi$","0",r"$\pi$",r"$2\pi$",r"$3\pi$",r"$4\pi$",r"$5\pi$"])
    Ikgcbar=fig.colorbar(Ikgsc, ax=axIkg)
    #Ikgcbar.set_ticks([-2*np.pi,-np.pi,0,np.pi,2*np.pi,3*np.pi,4*np.pi])
    #Ikgcbar.set_ticklabels([r"$-2\pi$",r"$-\pi$","0",r"$\pi$",r"$2\pi$",r"$3\pi$",r"$4\pi$"])
    Ik2cbar=fig.colorbar(Ik2sc, ax=axIk2)
    Itaucbar=fig.colorbar(Itausc, ax=axItau)

    #ax.set_xlabel(r"$\lambda\sigma_0/(k_B T)$")
    #ax.set_ylabel(r"$\tau_0\sigma_0$")
    fig.text(0.5,0.04,r"$B/(k_B T \sigma_0)$", ha="center")
    fig.text(0.04,0.5,r"$B'/(k_B T \sigma_0)$",va="center",rotation="vertical")
    #plt.tight_layout()
    #plt.show()
    plt.savefig(savefile, format="pdf",transparent=True)
    plt.close()