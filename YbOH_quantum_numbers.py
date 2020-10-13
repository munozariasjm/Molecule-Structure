import numpy as np

def q_numbers_174_bBJ(N_range,Lambda,S=1/2,I=1/2,all_M=True):
    Nmin,Nmax=N_range[0],N_range[-1]
    q_str = ['L','N','J','F','M']
    q_numbers = {}
    for q in q_str:
        q_numbers[q] = []
    if Lambda==0:
        for L in range(1):
            for N in np.arange(Nmin,Nmax+1,1):
                for J in np.arange(abs(N-S),abs(N+S)+1,1):
                    for F in np.arange(abs(J-I),abs(J+I)+1,1):
                        if all_M:
                            Mmin = -F
                        else:
                            Mmin = abs(F) % 1
                        for M in np.arange(Mmin,F+1,1):
                            values = [L,N,J,F,M]
                            for q,val in zip(q_str,values):
                                q_numbers[q].append(val+0)    #looks weird but adding 0 converts -0 to 0
    elif Lambda==1:
        if Nmin<abs(Lambda):
            print('Nmin must be >= L')
            Nmin=abs(Lambda)
        for N in np.arange(Nmin,Nmax+1,1):
            for J in np.arange(abs(N-S),abs(N+S)+1,1):
                for F in np.arange(abs(J-I),abs(J+I)+1,1):
                    if all_M:
                        Mmin = -F
                    else:
                        Mmin = abs(F) % 1
                    for M in np.arange(Mmin,F+1,1):
                        for L in range(-1,1+1,2):
                            values = [L,N,J,F,M]
                            for q,val in zip(q_str,values):
                                q_numbers[q].append(val+0)    #looks weird but adding 0 converts -0 to 0
    return q_numbers


def q_numbers_174_aBJ(N_range,Lambda=1,S=1/2,I=1/2,all_M=True,Omega_values=[1/2]):
    Nmin,Nmax=N_range[0],N_range[-1]
    Jmin = abs(Nmin-S)
    Jmax = abs(Nmax+S)
    q_str = ['L','Sigma','Omega','J','F','M']
    q_numbers = {}
    for q in q_str:
        q_numbers[q] = []
    if Lambda==0:
        for L in range(1):
            for J in np.arange(Jmin,Jmax+1,1):
                for F in np.arange(abs(J-I),abs(J+I)+1,1):
                    if all_M:
                        Mmin = -F
                    else:
                        Mmin = abs(F) % 1
                    for M in np.arange(Mmin,F+1,1):
                        for Sigma in np.arange(-abs(S),abs(S)+1,1):
                            Omega=L+Sigma
                            if abs(Omega) not in Omega_values:
                                continue
                            else:
                                values = [L,Sigma,Omega,J,F,M]
                            for q,val in zip(q_str,values):
                                q_numbers[q].append(val+0)    #looks weird but adding 0 converts -0 to 0
    elif Lambda==1:
        if Nmin<abs(Lambda):
            print('Nmin must be >= L')
            Nmin=abs(Lambda)
        Jmin = abs(Nmin-S)
        Jmax = abs(Nmax+S)
        for J in np.arange(Jmin,Jmax+1,1):
            for F in np.arange(abs(J-I),abs(J+I)+1,1):
                if all_M:
                    Mmin = -F
                else:
                    Mmin = abs(F) % 1
                for M in np.arange(Mmin,F+1,1):
                    for L in range(-Lambda,Lambda+1,2):
                        for Sigma in np.arange(-abs(S),abs(S)+1,1):
                            Omega=L+Sigma
                            if abs(Omega) not in Omega_values:
                                continue
                            else:
                                values = [L,Sigma,Omega,J,F,M]
                            for q,val in zip(q_str,values):
                                q_numbers[q].append(val+0)    #looks weird but adding 0 converts -0 to 0
    return q_numbers

def q_numbers_174_decoupled(N_range,Lambda=0, S=1/2, I=1/2,all_M = True):
    Nmin,Nmax=N_range[0],N_range[-1]
    q_str = ['L','N','M_N','M_S','M_I','M_F']
    q_numbers = {}
    for q in q_str:
        q_numbers[q] = []
    if Lambda==0:
        for L in range(1):
            for N in np.arange(Nmin,Nmax+1,1):
                for M_N in np.arange(-abs(N),abs(N)+1,1):
                    for M_S in np.arange(-abs(S),abs(S)+1,1):
                        for M_I in np.arange(-abs(I),abs(I)+1,1):
                            M_F = M_N + M_S + M_I
                            if (not all_M and M_F<0):
                                continue
                            else:
                                values = [L,N,M_N,M_S,M_I,M_F]
                                for q,val in zip(q_str,values):
                                    q_numbers[q].append(val+0)    #looks weird but adding 0 converts -0 to 0
    elif Lambda==1:
        if Nmin<abs(Lambda):
            print('Nmin must be >= L')
            Nmin=abs(Lambda)
        for N in np.arange(Nmin,Nmax+1,1):
            for M_N in np.arange(-abs(N),abs(N)+1,1):
                for M_S in np.arange(-abs(S),abs(S)+1,1):
                    for M_I in np.arange(-abs(I),abs(I)+1,1):
                        M_F = M_N + M_S + M_I
                        if (not all_M and M_F<0):
                            continue
                        else:
                            for L in range(-1,1+1,2):
                                values = [L,N,M_N,M_S,M_I,M_F]
                                for q,val in zip(q_str,values):
                                    q_numbers[q].append(val+0)    #looks weird but adding 0 converts -0 to 0
    return q_numbers



def q_numbers_173_aBJ(N_range,Lambda=1,S=1/2,I=5/2,iH=1/2,all_M=True,Omega_values=[1/2]):
    Nmin,Nmax=N_range[0],N_range[-1]
    q_str = ['L','Sigma','Omega','J','F1','F','M']
    q_numbers = {}
    for q in q_str:
        q_numbers[q] = []
    if Lambda==0:
        for L in range(1):
            for N in np.arange(Nmin,Nmax+1,1):
                for J in np.arange(abs(N-S),abs(N+S)+1,1):
                    for F1 in np.arange(abs(J-I),abs(J+I)+1,1):
                        for F in np.arange(abs(F1-iH),abs(F1+iH)+1,1):
                            if all_M:
                                Mmin = -F
                            else:
                                Mmin = abs(F) % 1
                            for M in np.arange(Mmin,F+1,1):
                                for Sigma in np.arange(-abs(S),abs(S)+1,1):
                                    Omega=L+Sigma
                                    if abs(Omega) not in Omega_values:
                                        continue
                                    else:
                                        values = [L,Sigma,Omega,J,F,M]
                                    for q,val in zip(q_str,values):
                                        q_numbers[q].append(val+0)    #looks weird but adding 0 converts -0 to 0
    elif Lambda==1:
        if Nmin<abs(Lambda):
            print('Nmin must be >= L')
            Nmin=abs(Lambda)
        for N in np.arange(Nmin,Nmax+1,1):
            for J in np.arange(abs(N-S),abs(N+S)+1,1):
                for F1 in np.arange(abs(J-I),abs(J+I)+1,1):
                    for F in np.arange(abs(F1-iH),abs(F1+iH)+1,1):
                        if all_M:
                            Mmin = -F
                        else:
                            Mmin = abs(F) % 1
                        for M in np.arange(Mmin,F+1,1):
                            for L in range(-Lambda,Lambda+1,2):
                                for Sigma in np.arange(-abs(S),abs(S)+1,1):
                                    Omega=L+Sigma
                                    if abs(Omega) not in Omega_values:
                                        continue
                                    else:
                                        values = [L,Sigma,Omega,J,F,M]
                                    for q,val in zip(q_str,values):
                                        q_numbers[q].append(val+0)    #looks weird but adding 0 converts -0 to 0
    return q_numbers

def q_numbers_173_bBS(N_range,Lambda,S=1/2,I=5/2,iH=1/2,all_M=True):
    Nmin,Nmax=N_range[0],N_range[-1]
    q_str = ['L','N','G','F1','F','M']
    q_numbers = {}
    for q in q_str:
        q_numbers[q] = []
    if Lambda==0:
        for L in range(1):
            for N in np.arange(Nmin,Nmax+1,1):
                for G in np.arange(abs(I-S),abs(I+S)+1,1):
                    for F1 in np.arange(abs(G-N),abs(G+N)+1,1):
                        for F in np.arange(abs(F1-iH),abs(F1+iH)+1,1):
                            if all_M:
                                Mmin = -F
                            else:
                                Mmin = abs(F) % 1
                            for M in np.arange(Mmin,F+1,1):
                                values = [L,N,G,F1,F,M]
                                for q,val in zip(q_str,values):
                                    q_numbers[q].append(val+0)    #looks weird but adding 0 converts -0 to 0
    elif Lambda==1:
        if Nmin<abs(Lambda):
            print('Nmin must be >= L')
            Nmin=abs(Lambda)
        for N in np.arange(Nmin,Nmax+1,1):
            for G in np.arange(abs(I-S),abs(I+S)+1,1):
                for F1 in np.arange(abs(G-N),abs(G+N)+1,1):
                    for F in np.arange(abs(F1-iH),abs(F1+iH)+1,1):
                        if all_M:
                            Mmin = -F
                        else:
                            Mmin = abs(F) % 1
                        for M in np.arange(Mmin,F+1,1):
                            for L in range(-1,1+1,2):
                                values = [L,N,G,F1,F,M]
                                for q,val in zip(q_str,values):
                                    q_numbers[q].append(val+0)    #looks weird but adding 0 converts -0 to 0
    return q_numbers


def q_numbers_173_decoupled(N_range,Lambda=0, S=1/2, I=5/2, iH=1/2, all_M = True):
    Nmin,Nmax=N_range[0],N_range[-1]
    q_str = ['L','N','M_N','M_S','M_I','M_iH','M_F']
    q_numbers = {}
    for q in q_str:
        q_numbers[q] = []
    if Lambda==0:
        for L in range(1):
            for N in np.arange(Nmin,Nmax+1,1):
                for M_N in np.arange(-abs(N),abs(N)+1,1):
                    for M_S in np.arange(-abs(S),abs(S)+1,1):
                        for M_I in np.arange(-abs(I),abs(I)+1,1):
                            for M_iH in np.arange(-abs(iH),abs(iH)+1,1):
                                M_F = M_N + M_S + M_I+ M_iH
                                if (not all_M and M_F<0):
                                    continue
                                else:
                                    values = [L,N,M_N,M_S,M_I,M_iH,M_F]
                                    for q,val in zip(q_str,values):
                                        q_numbers[q].append(val+0)    #looks weird but adding 0 converts -0 to 0
    elif Lambda==1:
        if Nmin<abs(Lambda):
            print('Nmin must be >= L')
            Nmin=abs(Lambda)
        for N in np.arange(Nmin,Nmax+1,1):
            for M_N in np.arange(-abs(N),abs(N)+1,1):
                for M_S in np.arange(-abs(S),abs(S)+1,1):
                    for M_I in np.arange(-abs(I),abs(I)+1,1):
                        for M_iH in np.arange(-abs(iH),abs(iH)+1,1):
                            M_F = M_N + M_S + M_I+ M_iH
                            if (not all_M and M_F<0):
                                continue
                            else:
                                for L in range(-1,1+1,2):
                                    values = [L,N,M_N,M_S,M_I,M_iH,M_F]
                                    for q,val in zip(q_str,values):
                                        q_numbers[q].append(val+0)    #looks weird but adding 0 converts -0 to 0
    return q_numbers
