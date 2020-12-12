import numpy as np
import matplotlib.pyplot as plt
from plot import *
from autocorrelation import *


def read_multi_file(foldername, kar, lam, m):
    En = []
    En_e = []
    lam_L_e = []
    for i in range(m + 1)[1:]:
        filename = foldername + \
            "/E_v_MC_kar" + str(kar) + "_lam" + \
            str(lam) + "_%de5.txt" % m
        En_b, En_e_b, lam_L_e_b = np.loadtxt(
            filename, skiprows=1, delimiter=",", unpack=True)
        En = np.concatenate([En, En_b])
        En_e = np.concatenate([En_e, En_e_b])
        lam_L_e = np.concatenate([lam_L_e, lam_L_e_b])
    return (En, En_e, lam_L_e)

def O_stat_ana(foldername, N_, Ne_,L_, kar_,C0_, karg_,lam_, B_, Bc_, tau0_ ,c,mode):
    E_ave, E_tau, E_err = [], [], []
    Les_ave, Les_tau, Les_err = [[] for i in range(Ne_)], [[] for i in range(Ne_)], [[] for i in range(Ne_)]
    IdA_ave, IdA_tau, IdA_err = [], [], []
    I2H_ave, I2H_tau, I2H_err = [], [], []
    I2H2_ave, I2H2_tau, I2H2_err = [], [], []
    Ikg_ave, Ikg_tau, Ikg_err = [], [], []
    Ik2_ave, Ik2_tau, Ik2_err = [], [], []
    Itau_ave, Itau_tau, Itau_err = [], [], []
    #Itau2_ave, Itau2_tau, Itau2_err = [], [], []
    Tau_ave, Tau_tau, Tau_err = [], [], []  # <Itau/Le>
    Ieig0_ave,Ieig0_tau,Ieig0_err, = [],[],[] #smallest eigenvalue
    Ieig1_ave,Ieig1_tau,Ieig1_err, = [],[],[] #eigenvalue
    Ieig2_ave,Ieig2_tau,Ieig2_err, = [],[],[] #largest eigenvalue
    if(Ne_==2):
        Ledif_ave,Ledif_tau,Ledif_err=[],[],[]
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
            N,Ne, L, kar, C0, karg, lam, B, Bc, tau0 = N_, Ne_,L_, kar_,C0_, karg_, lam_[i], B_, Bc_, tau0_
        elif (mode == "tau0"):
            N,Ne, L, kar, C0, karg, lam, B, Bc, tau0 = N_, Ne_,L_, kar_, C0_,karg_,lam_,  B_, Bc_, tau0_[
                i]
        elif (mode == "N"):
            N,Ne, L, kar, C0, karg, lam, B, Bc, tau0 = N_[i], Ne_, L_, kar_,C0_, karg_,lam_, B_, Bc_, tau0_
        elif (mode == "L"):
            N,Ne, L, kar,C0,  karg,lam, B, Bc, tau0 = N_, Ne_,L_[i], kar_,C0_, karg_,lam_,  B_, Bc_, tau0_
        elif (mode == "B"):
            N,Ne, L, kar, C0, karg,lam, B, Bc, tau0 = N_, Ne_,L_, kar_, C0_, karg_,lam_,  B_[i], Bc_, tau0_
        elif (mode == "Bc"):
            N,Ne, L, kar, C0, karg, lam, B, Bc, tau0 = N_,Ne_, L_, kar_, C0_, karg_,lam_,  B_, Bc_[i], tau0_
        print("dealing with ", N,Ne, L, kar, C0, karg, lam, B, Bc, tau0)

        f2rtail = "MC_N%.0f_Ne%.0f_L%.0f_kar%.0f_C0%.1f_karg%.1f_lam%.1f_B%.1f_Bc%.0f_tau0%.2f.txt" % (N, Ne, L, kar, C0, karg, lam, B, Bc, tau0)
        file2read = foldername + "/O_"+f2rtail
        data = np.loadtxt(file2read, skiprows=13, delimiter=",", unpack=True)
        E = data[0]
        Les = data[1:1+Ne]
        IdA,I2H,I2H2, Ikg, Ik2, Itau = data[1+Ne:]
        file2readIij = foldername + "/Iij_"+f2rtail
        Ieig0,Ieig1,Ieig2 = [[0],[0],[0]]
        Ieig0,Ieig1,Ieig2 = data_eigen_Iij(file2readIij)


        # Ne2 case, need Ledif for additional info
        if(Ne==2):
            Ledif = np.abs(Les[0]-Les[1])
            Ledif_ave.append(np.average(Ledif))
            rho, cov0 = autocorrelation_function_fft(Ledif)
            tau, tau_err = tau_int_cal_rho(rho,c)
            Ledif_tau.append(tau)
            Ledif_err.append(np.sqrt(2 * tau / len(Ledif) * cov0))

        # E
        E_ave.append(np.average(E))
        rho, cov0 = autocorrelation_function_fft(E)
        #print("cov0-var(E)",cov0-np.var(E))
        tau, tau_err = tau_int_cal_rho(rho,c)
        # autocorrelation_plot(rho, tau, file2read[:-4] + "_autoE.pdf")
        E_tau.append(tau)
        E_err.append(np.sqrt(2 * tau / len(E) * cov0))

        # Le
        for e in range(Ne):
            Les_ave[e].append(np.average(Les[e]))
            rho, cov0 = autocorrelation_function_fft(Les[e])
            tau, tau_err = tau_int_cal_rho(rho,c)
            # autocorrelation_plot(rho, tau, file2read[:-4] + "_autoLe.pdf")
            Les_tau[e].append(tau)
            Les_err[e].append(np.sqrt(2 * tau / len(Les[e]) * cov0))

        # IdA
        IdA_ave.append(np.average(IdA))
        rho, cov0 = autocorrelation_function_fft(IdA)
        tau, tau_err = tau_int_cal_rho(rho,c)
        IdA_tau.append(tau)
        IdA_err.append(np.sqrt(2 * tau / len(IdA) * cov0))

        # I2H
        I2H_ave.append(np.average(I2H))
        rho, cov0 = autocorrelation_function_fft(I2H)
        tau, tau_err = tau_int_cal_rho(rho,c)
        I2H_tau.append(tau)
        I2H_err.append(np.sqrt(2 * tau / len(I2H) * cov0))

        # I2H2
        I2H2_ave.append(np.average(I2H2))
        rho, cov0 = autocorrelation_function_fft(I2H2)
        tau, tau_err = tau_int_cal_rho(rho,c)
        # autocorrelation_plot(rho, tau, file2read[:-4] + "_autoI2H2.pdf")
        I2H2_tau.append(tau)
        I2H2_err.append(np.sqrt(2 * tau / len(I2H2) * cov0))

        # Ikg
        Ikg_ave.append(np.average(Ikg))
        rho, cov0 = autocorrelation_function_fft(Ikg)
        tau, tau_err = tau_int_cal_rho(rho,c)
        # autocorrelation_plot(rho, tau, file2read[:-4] + "_autoIkg.pdf")
        Ikg_tau.append(tau)
        Ikg_err.append(np.sqrt(2 * tau / len(Ikg) * cov0))

        # Ik2
        Ik2_ave.append(np.average(Ik2))
        rho, cov0 = autocorrelation_function_fft(Ik2)
        tau, tau_err = tau_int_cal_rho(rho,c)
        # autocorrelation_plot(rho, tau, file2read[:-4] + "_autoIk2.pdf")
        Ik2_tau.append(tau)
        Ik2_err.append(np.sqrt(2 * tau / len(Ik2) * cov0))

        # Itau
        Itau_ave.append(np.average(Itau))
        rho, cov0 = autocorrelation_function_fft(Itau)
        tau, tau_err = tau_int_cal_rho(rho,c)
        # autocorrelation_plot(rho, tau, file2read[:-4] + "_autoItau.pdf")
        Itau_tau.append(tau)
        Itau_err.append(np.sqrt(2 * tau / len(Itau) * cov0))

        # Ieig
        Ieig0_ave.append(np.average(Ieig0))
        rho, cov0 = autocorrelation_function_fft(Ieig0)
        tau, tau_err = tau_int_cal_rho(rho,c)
        Ieig0_tau.append(tau)
        Ieig0_err.append(np.sqrt(2*tau/len(Ieig0)*cov0))

        Ieig1_ave.append(np.average(Ieig1))
        rho, cov0 = autocorrelation_function_fft(Ieig1)
        tau, tau_err = tau_int_cal_rho(rho,c)
        Ieig1_tau.append(tau)
        Ieig1_err.append(np.sqrt(2*tau/len(Ieig1)*cov0))

        Ieig2_ave.append(np.average(Ieig2))
        rho, cov0 = autocorrelation_function_fft(Ieig2)
        tau, tau_err = tau_int_cal_rho(rho,c)
        Ieig2_tau.append(tau)
        Ieig2_err.append(np.sqrt(2*tau/len(Ieig2)*cov0))


    # only changed "lam" and "B" mode here, others waiting for further decision
    if(mode == "lam"):
        f2stail = "MC_N%.0f_Ne%.0f_L%.0f_kar%.0f_C0%.1f_karg%.1f_lams_B%.1f_Bc%.0f_tau0%.2f_ana.txt" % (N, Ne,L, kar, C0, karg, B, Bc, tau0)
    elif(mode == "tau0"):
        f2stail = "MC_N%.0f_Ne%.0f_L%.0f_kar%.0f_lam%.1f_karg%.1f_B%.1f_Bc%.0f_tau0s_ana.txt" % (N, Ne,L, kar, lam, karg, B, Bc)
    elif(mode == "N"):
        f2stail = "MC_Ns_Ne%.0f_L%.0f_kar%.0f_lam%.1f_karg%.1f_B%.1f_Bc%.0f_tau0%.2f_ana.txt" % (L, kar, lam, karg, B, Bc, tau0)
    elif(mode == "L"):
        f2stail = "MC_N%.0f_Ne%.0f_Ls_kar%.0f_lam%.1f_karg%.1f_B%.1f_Bc%.0f_tau0%.2f_ana.txt" % (N, Ne,kar, lam, karg, B, Bc, tau0)
    elif(mode == "B"):
        f2stail = "MC_N%.0f_Ne%.0f_L%.0f_kar%.0f_C0%.1f_karg%.1f_lam%.1f_Bs_Bc%.0f_tau0%.2f_ana.txt" % (N, Ne,L, kar, C0, karg, lam, Bc, tau0)
    elif(mode == "Bc"):
        f2stail = "MC_N%.0f_Ne%.0f_L%.0f_kar%.0f_lam%.1f_karg%.1f_B%.1f_Bcs_tau0%.2f_ana.txt" % (N, Ne,L, kar, lam, karg, B, tau0)

    savefile = foldername + "/O_" + f2stail
    savefileIij = foldername + "/Iij_" + f2stail
    with open(savefile, "w") as f:
        f.write(mode+",E_ave,E_tau,E_err")
        for e in range(Ne):
            f.write(",Les_ave[%d],Les_tau[%d],Les_err[%d]"%(e,e,e))
        f.write(",IdA_ave,IdA_tau,IdA_err,I2H_ave,I2H_tau,I2H_err,I2H2_ave,I2H2_tau,I2H2_err,Ikg_ave,Ikg_tau,Ikg_err,Ik2_ave,Ik2_tau,Ik2_err,Itau_ave,Itau_tau,Itau_err")
        if(Ne==2):
            f.write(",Ledif_ave,Ledif_tau,Ledif_err")
        f.write("\n")
        for i in range(len(cpar)):
            f.write("%f,%f,%f,%f" % (cpar[i], E_ave[i], E_tau[i], E_err[i]))
            for e in range(Ne):
                f.write(",%f,%f,%f"%(Les_ave[e][i],Les_tau[e][i], Les_err[e][i]))
            f.write(",%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f"%(IdA_ave[i], IdA_tau[i], IdA_err[i],I2H_ave[i], I2H_tau[i], I2H_err[i],I2H2_ave[i], I2H2_tau[i], I2H2_err[i], Ikg_ave[i], Ikg_tau[i], Ikg_err[i], Ik2_ave[i], Ik2_tau[i], Ik2_err[i], Itau_ave[i], Itau_tau[i], Itau_err[i]))
            if(Ne==2):
                f.write(",%f,%f,%f"%(Ledif_ave[i], Ledif_tau[i], Ledif_err[i]))
            f.write("\n")


    with open(savefileIij,"w") as fIij:
        fIij.write(mode+",Ieig0_ave,Ieig0_tau,Ieig0_err,Ieig1_ave,Ieig1_tau,Ieig1_err,Ieig2_ave,Ieig2_tau,Ieig2_err\n")
        for i in range(len(cpar)):
            fIij.write("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n"%(cpar[i],Ieig0_ave[i],Ieig0_tau[i],Ieig0_err[i],Ieig1_ave[i],Ieig1_tau[i],Ieig1_err[i],Ieig2_ave[i],Ieig2_tau[i],Ieig2_err[i]))


def data_eigen_Iij(filename):
    Iijs = np.loadtxt(filename,delimiter=",",skiprows=1)
    Ieigs = []
    for Iij in Iijs:
        #print("Iij",Iij)
        #print(np.reshape(Iij,(3,3)))
        w,v = np.linalg.eig(np.reshape(Iij,(3,3)))
        Ieigs.append(np.sort(w))
    Ieigs = np.transpose(Ieigs)
    return Ieigs


def O_stat_ana_Ls(foldername, L, kar, lam):
    En_ave, En_std = [], []
    L_ave, L_std = [], []
    for i in range(len(L)):
        config_plot_xyz(foldername + "/State_L" +
                        str(L[i]) + "_kar" + str(kar) + "_lam" + str(lam) + ".txt")
        print("dealing with ", L[i], kar, lam)
        file2read = foldername + "/O_MC_L" + \
            str(L[i]) + "_kar" + str(kar) + "_lam" + str(lam) + ".txt"
        En, L_e = np.loadtxt(
            file2read, skiprows=7, delimiter=",", unpack=True)
        i2H = (En - lam * L_e) / kar
        O_kar_lam_MCstep_plot(L_e, i2H, En, file2read[:-4] + ".pdf")
        L_ave.append(np.average(L_e))
        L_std.append(np.std(L_e))
        En_ave.append(np.average(En))
        En_std.append(np.std(En))
    savefile = foldername + "/O_MC_kar" + \
        str(kar) + "_lam" + str(lam) + "ana.txt"
    with open(savefile, "w") as f:
        f.write("L,En_ave,En_std,L_ave,L_std\n")
        for i in range(len(L)):
            f.write("%f,%f,%f,%f,%f\n" % (
                L[i], En_ave[i], En_std[i], L_ave[i], L_std[i]))
