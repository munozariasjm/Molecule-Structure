import numpy as np


def q_numbers_174_bBJ(N_range,Lambda,S=1/2,I_list=[0,1/2],M_values='all'):
    I=I_list[-1]
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
                        if M_values=='none':
                            M=abs(F) % 1
                            values = [L,N,J,F,M]
                            for q,val in zip(q_str,values):
                                q_numbers[q].append(val+0)    #looks weird but adding 0 converts -0 to 0
                        else:
                            if M_values=='all':
                                Mmin = -F
                            elif M_values=='pos':
                                Mmin = abs(F) % 1
                            for M in np.arange(Mmin,F+1,1):
                                values = [L,N,J,F,M]
                                for q,val in zip(q_str,values):
                                    q_numbers[q].append(val+0)    #looks weird but adding 0 converts -0 to 0
    elif Lambda>0:
        if Nmin<abs(Lambda):
            print('Nmin must be >= L')
            Nmin=abs(Lambda)
        for N in np.arange(Nmin,Nmax+1,1):
            for J in np.arange(abs(N-S),abs(N+S)+1,1):
                for F in np.arange(abs(J-I),abs(J+I)+1,1):
                    if M_values=='none':
                        for L in range(-Lambda,Lambda+1,2):
                            M=abs(F) % 1
                            values = [L,N,J,F,M]
                            for q,val in zip(q_str,values):
                                q_numbers[q].append(val+0)    #looks weird but adding 0 converts -0 to 0
                    else:
                        if M_values=='all':
                            Mmin = -F
                        elif M_values=='pos':
                            Mmin = abs(F) % 1
                        for M in np.arange(Mmin,F+1,1):
                            for L in range(-Lambda,Lambda+1,2):
                                values = [L,N,J,F,M]
                                for q,val in zip(q_str,values):
                                    q_numbers[q].append(val+0)    #looks weird but adding 0 converts -0 to 0
    return q_numbers


def q_numbers_174_aBJ(N_range,Lambda=1,S=1/2,I_list=[0,1/2],M_values='all',Omega_values=[1/2]):
    I=I_list[-1]
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
                    if M_values=='none':
                        M=abs(F) % 1
                        for Sigma in np.arange(-abs(S),abs(S)+1,1):
                            Omega=L+Sigma
                            if abs(Omega) not in Omega_values:
                                continue
                            else:
                                values = [L,Sigma,Omega,J,F,M]
                            for q,val in zip(q_str,values):
                                q_numbers[q].append(val+0)    #looks weird but adding 0 converts -0 to 0
                    else:
                        if M_values=='all':
                            Mmin = -F
                        elif M_values=='pos':
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
    elif Lambda>0:
        if Nmin<abs(Lambda):
            print('Nmin must be >= L')
            Nmin=abs(Lambda)
        Jmin = abs(Nmin-S)
        Jmax = abs(Nmax+S)
        for J in np.arange(Jmin,Jmax+1,1):
            for F in np.arange(abs(J-I),abs(J+I)+1,1):
                if M_values=='none':
                    M=abs(F) % 1
                    for L in range(-Lambda,Lambda+1,2):
                        for Sigma in np.arange(-abs(S),abs(S)+1,1):
                            Omega=L+Sigma
                            if abs(Omega) not in Omega_values:
                                continue
                            else:
                                values = [L,Sigma,Omega,J,F,M]
                            for q,val in zip(q_str,values):
                                q_numbers[q].append(val+0)    #looks weird but adding 0 converts -0 to 0
                else:
                    if M_values=='all':
                        Mmin = -F
                    elif M_values=='pos':
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

def q_numbers_174_decoupled(N_range,Lambda=0, S=1/2, I_list=[0,1/2],M_values='all'):
    I = I_list[-1]
    if M_values != 'all':
        print('Cannot construct decoupled basis without M values')
        return {}
    else:
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
                                if (M_values!='all' and M_F<0):
                                    continue
                                else:
                                    values = [L,N,M_N,M_S,M_I,M_F]
                                    for q,val in zip(q_str,values):
                                        q_numbers[q].append(val+0)    #looks weird but adding 0 converts -0 to 0
        elif Lambda>1:
            if Nmin<abs(Lambda):
                print('Nmin must be >= L')
                Nmin=abs(Lambda)
            for N in np.arange(Nmin,Nmax+1,1):
                for M_N in np.arange(-abs(N),abs(N)+1,1):
                    for M_S in np.arange(-abs(S),abs(S)+1,1):
                        for M_I in np.arange(-abs(I),abs(I)+1,1):
                            M_F = M_N + M_S + M_I
                            if (M_values!='all' and M_F<0):
                                continue
                            else:
                                for L in range(-Lambda,Lambda+1,2):
                                    values = [L,N,M_N,M_S,M_I,M_F]
                                    for q,val in zip(q_str,values):
                                        q_numbers[q].append(val+0)    #looks weird but adding 0 converts -0 to 0
    return q_numbers



def q_numbers_173_aBJ(N_range,Lambda=1,S=1/2,I_list=[5/2,1/2],M_values='all',Omega_values=[1/2]):
    I = I_list[0]
    iH = I_list[-1]
    Nmin,Nmax=N_range[0],N_range[-1]
    Jmin = abs(Nmin-S)
    Jmax = abs(Nmax+S)
    q_str = ['L','Sigma','Omega','J','F1','F','M']
    q_numbers = {}
    for q in q_str:
        q_numbers[q] = []
    if Lambda==0:
        for L in range(1):
            for J in np.arange(Jmin,Jmax+1,1):
                for F1 in np.arange(abs(J-I),abs(J+I)+1,1):
                    for F in np.arange(abs(J-iH),abs(J+iH)+1,1):
                        if M_values=='none':
                            for Sigma in np.arange(-abs(S),abs(S)+1,1):
                                Omega=L+Sigma
                                if abs(Omega) not in Omega_values:
                                    continue
                                else:
                                    values = [L,Sigma,Omega,J,F1,F]
                                for q,val in zip(q_str,values):
                                    q_numbers[q].append(val+0)    #looks weird but adding 0 converts -0 to 0
                        else:
                            if M_values=='all':
                                Mmin = -F
                            elif M_values=='pos':
                                Mmin = abs(F) % 1
                            for M in np.arange(Mmin,F+1,1):
                                for Sigma in np.arange(-abs(S),abs(S)+1,1):
                                    Omega=L+Sigma
                                    if abs(Omega) not in Omega_values:
                                        continue
                                    else:
                                        values = [L,Sigma,Omega,J,F1,F,M]
                                    for q,val in zip(q_str,values):
                                        q_numbers[q].append(val+0)    #looks weird but adding 0 converts -0 to 0
    elif Lambda>0:
        if Nmin<abs(Lambda):
            print('Nmin must be >= L')
            Nmin=abs(Lambda)
        Jmin = abs(Nmin-S)
        Jmax = abs(Nmax+S)
        for N in np.arange(Nmin,Nmax+1,1):
            for J in np.arange(abs(N-S),abs(N+S)+1,1):
                for F1 in np.arange(abs(J-I),abs(J+I)+1,1):
                    for F in np.arange(abs(F1-iH),abs(F1+iH)+1,1):
                        if M_values=='none':
                            for Sigma in np.arange(-abs(S),abs(S)+1,1):
                                for L in range(-Lambda,Lambda+1,2):
                                    Omega=L+Sigma
                                    if abs(Omega) not in Omega_values:
                                        continue
                                    else:
                                        values = [L,Sigma,Omega,J,F1,F]
                                    for q,val in zip(q_str,values):
                                        q_numbers[q].append(val+0)    #looks weird but adding 0 converts -0 to 0
                        else:
                            if M_values=='all':
                                Mmin = -F
                            elif M_values=='pos':
                                Mmin = abs(F) % 1
                            for M in np.arange(Mmin,F+1,1):
                                for L in range(-Lambda,Lambda+1,2):
                                    for Sigma in np.arange(-abs(S),abs(S)+1,1):
                                        Omega=L+Sigma
                                        if abs(Omega) not in Omega_values:
                                            continue
                                        else:
                                            values = [L,Sigma,Omega,J,F1,F,M]
                                        for q,val in zip(q_str,values):
                                            q_numbers[q].append(val+0)    #looks weird but adding 0 converts -0 to 0

    return q_numbers

def q_numbers_173_bBS(N_range,Lambda,S=1/2,I_list=[5/2,1/2],M_values='all'):
    I=I_list[0]
    iH=I_list[-1]
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
                            if M_values=='none':
                                M= abs(F) % 1
                                values = [L,N,G,F1,F,M]
                                for q,val in zip(q_str,values):
                                    q_numbers[q].append(val+0)    #looks weird but adding 0 converts -0 to 0
                            else:
                                if M_values=='all':
                                    Mmin = -F
                                elif M_values=='pos':
                                    Mmin = abs(F) % 1
                                for M in np.arange(Mmin,F+1,1):
                                    values = [L,N,G,F1,F,M]
                                    for q,val in zip(q_str,values):
                                        q_numbers[q].append(val+0)    #looks weird but adding 0 converts -0 to 0
    elif Lambda>1:
        if Nmin<abs(Lambda):
            print('Nmin must be >= L')
            Nmin=abs(Lambda)
        for N in np.arange(Nmin,Nmax+1,1):
            for G in np.arange(abs(I-S),abs(I+S)+1,1):
                for F1 in np.arange(abs(G-N),abs(G+N)+1,1):
                    for F in np.arange(abs(F1-iH),abs(F1+iH)+1,1):
                        if M_values=='none':
                            for L in range(-Lambda,Lambda+1,2):
                                M=abs(F) % 1
                                values = [L,N,G,F1,F,M]
                                for q,val in zip(q_str,values):
                                    q_numbers[q].append(val+0)    #looks weird but adding 0 converts -0 to 0
                        else:
                            if M_values=='all':
                                Mmin = -F
                            elif M_values=='pos':
                                Mmin = abs(F) % 1
                            for M in np.arange(Mmin,F+1,1):
                                for L in range(-Lambda,Lambda+1,2):
                                    values = [L,N,G,F1,F,M]
                                    for q,val in zip(q_str,values):
                                        q_numbers[q].append(val+0)    #looks weird but adding 0 converts -0 to 0
    return q_numbers


def q_numbers_173_decoupled(N_range,Lambda=0, S=1/2, I_list=[5/2,1/2],M_values='all'):
    I = I_list[0]
    iH = I_list[-1]
    if M_values != 'all':
        print('Cannot construct decoupled basis without M values')
        return {}
    else:
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
                                    if (M_values!='all' and M_F<0):
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
                                if (M_values!='all' and M_F<0):
                                    continue
                                else:
                                    for L in range(-1,1+1,2):
                                        values = [L,N,M_N,M_S,M_I,M_iH,M_F]
                                        for q,val in zip(q_str,values):
                                            q_numbers[q].append(val+0)    #looks weird but adding 0 converts -0 to 0
        return q_numbers
