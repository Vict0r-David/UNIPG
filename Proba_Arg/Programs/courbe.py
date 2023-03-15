import numpy as np
import matplotlib.pyplot as plt

def c1():
    x1 = [10, 15, 16, 17, 18, 19, 20, 21]
    x2 = [10, 15, 16, 17, 18, 19, 20, 21, 30, 40, 50, 60, 70, 80, 90]
    constellation = [0.06, 2.2, 4.5, 11, 20.5, 44.5, 78, 176]
    avg_MCN = [0.0003, 0.0004, 0.0004, 0.0004, 0.0004, 0.0004, 0.0004, 0.0008, 0.002, 0.025, 0.004, 0.035, 0.03, 0.3, 1]
    max_MCN = [0.05, 0.06, 0.05, 0.07, 0.05, 0.02, 0.05, 0.05, 0.06, 0.05, 0.25, 10, 3, 55, 198]

    plt.plot(x1, constellation,"b-+", label="Constellation")
    plt.plot(x2, avg_MCN,"g->", label="avg MCN")
    plt.plot(x2, max_MCN,"r-o", label="max MCN")

    plt.legend()
    plt.xlabel("Number of attacks")
    plt.ylabel("Time (sec)")

    plt.show()


def c2():
    x = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
    max_MCN = [0.15, 0.44, 1, 1.65, 2, 4, 5.1, 7.6, 10, 12.6]
    max_SCN = [0.03, 0.03, 0.03, 0.032, 0.0325, 0.0385, 0.04, 0.0425, 0.061, 0.071]

    diff_MCN = [0.43154666, 0.022956907, 0.026399593, 0.167250288, 0.043337491, 0.470326027, 0.122668016, 0.373260442, 0.53676511, 0.2612704]
    diff_SCN = [1.805253425, 0.022956907, 1.566724965, 2.071470837, 2.696645459, 0.869836158, 6.206188313, 13.17226876, 6.845911462, 11.22025839]

    #fig, host = plt.subplots(figsize=(8,5), layout='constrained')
    #ax2 = host.twinx()
    
    #host.set_xlim(0, 2)
    #host.set_ylim(0, 2)
    #ax2.set_ylim(0, 1)
        
    #host.set_xlabel("Number of attacks")
    #host.set_ylabel("Time (sec)")
    #ax2.set_ylabel("% Error")

    #color1, color2, color3 = plt.cm.viridis([0, .5, .9])

    #p1 = host.plot(x, diff_SCN, "r:^", label="% error M-C / SCN")
    #p2 = ax2.plot( x, diff_MCN, "b:>", label="% error M-C / MCN")

    #host.legend(handles=p1+p2, loc='best')


    #plt.plot(x1, constellation,"b-+", label="Constellation")
    plt.plot(x, max_SCN,"b-*", label="max time SCN")
    plt.plot(x, max_MCN,"r-o", label="max time MCN")

    plt.plot(x, diff_SCN,"b:^", label="% error M-C / SCN")
    plt.plot( x, diff_MCN,"r:>", label="% error M-C / MCN")

    plt.legend()
    plt.xlabel("Number of attacks")
    plt.ylabel("Time (sec) and % Error")
    

    plt.show()

#c1()
c2()