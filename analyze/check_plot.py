import numpy as np
import matplotlib.pyplot as plt


def curvature_compare():
    print("checking curvature")

    theta = np.pi * np.array([1 / 8, 1 / 6, 1 / 4, 1 / 3])
    thetap = np.pi * np.linspace(0, 0.5, 100)
    Eh = np.array([1.80031, 3.16238, 6.86592, 11.548])
    Eg = np.array([0.450093, 0.790644, 1.71672, 2.8877])
    plt.figure()
    plt.plot(theta, Eh, marker="s", linestyle="None",
             color="tomato", label="bending energy")
    plt.plot(thetap, 4 * 2 * np.pi * (1 - np.cos(thetap)), marker="None", linestyle="-",
             color="tomato", label=r"$8\pi(1-\cos(\theta))$")
    plt.plot(theta, Eg, marker="d", linestyle="None",
             color="royalblue", label="gaussian energy")
    plt.plot(thetap, 2 * np.pi * (1 - np.cos(thetap)), marker="None", linestyle="-",
             color="royalblue", label=r"$2\pi(1-\cos(\theta))$")
    plt.xlabel(r"$\theta$")
    plt.ylabel("energy")
    plt.legend()
    # plt.show()
    plt.savefig("cap_energy.pdf", dpi=300, format="pdf")
    plt.close()
    print(theta)
