 # Test:
    if name=='test':
        # Artificial data is constructed to test the softwares!
        t_int = [0, 200]                                       # [days]
        f_int = [0, 100]                                       # [c/d]
        t     = linspace(t_int[0], t_int[1], 1e4)
        # Stellar signal:
        f      = [0.005, 0.035, 0.003, 0.08]                   # Frequency [c/d] 
        A      = [700, 500, 500, 300]                          # Amplitude [Arbt. unit]
        Phi    = [0, pi/2, pi/3, pi/10]                        # Phase
        S_star = sum([A[i]*sin(2*pi*f[i]*t + Phi[i]) for i in range(len(f))], 0)
        # Slow trend: 
        S_slow = array(S_star[0]-range(len(S_star)))*0.2       # Linear function is added
        S_sum0 = S_star + S_slow
        # Signal with noise:
        S_noise = np.random.normal(t_int[0], t_int[1], 1e4)/5  # Random noise signal
        S_sum1  = S_sum0 + S_noise  
        # plt.plot(t, S_star, 'k-')
        # plt.plot(t, S_sum1, 'r-')
        # plt.plot(t, S_sum0, 'b-')
        # plt.show(); sys.exit()
        Pf, _, _, _, = power(vstack([t, S_sum1]).T, [0, 0.5], 1e-4)
        plot_multi('power', 'f', 'P', Pf)
        sys.exit()
        #---- Corrections:
        data = vstack([t, S_sum1]).T
        # data = locate(data, None, 3, 1)                 # Corrects for bad data
        data = jumps(data, 300, 1)                           # Corrects for jumps
        data = stellar_noise(data, [0.02, 1], clean, 2, 1) # Cleans for stellar signals.
        data = slowtrend(data, 10, 10, 1)                    # Corrects for slow trends.
        # # Exoplanet signal:
        # Psi    = [3]                               # Phase [hours]
        # P      = [80]                              # 2 transits [hours]
        # T      = 2.5
        # data_exo = loadtxt(os.path.join('/home/nicholas/Dropbox/Uni/Timeseries/Data/Exo_bach')) 
        # t_exo  = data_exo[:,0]
        # S_exo  = -data_exo[:,2] 
        # t_exo = t_exo + P*t_exo + 2*P*t_exo
        # # plt.plot(t_exo, S_exo); plt.show()
        # print len(t), len(t_exo), len(S_star), len(S_exo)
        # sys.exit()
        # S_sum0 = S_star + S_exo
        # plt.plot(t, S_sum0)
        # Slow trend:
   
