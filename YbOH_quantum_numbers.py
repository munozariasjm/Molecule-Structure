import numpy as np

def count_momenta(J_list, M_values):
    if M_values == 'all':
         num = 1
         for J in J_list:
             num*=(2*J+1)
    else:
        num = recursive_count(J_list,M_values,1)
    return num

def add_J(J1,J2):
    nsum = J1+J2 - abs(J1-J2) + 1
    Jsum = np.arange(abs(J1-J2),J1+J2+1, 1)
    return nsum,Jsum

def add_J_list(J_list):
    for (i,J) in enumerate(J_list):


def recursive_count(J_list,M_values,i,J = None):
    J0 = {True:J_list[0], False:J}[J is None]
    Ji = J_list[i]
    n = len(J_list) - i
    sum = 0
    if if n>1:
        for Jf in np.arange(abs(J0-Ji),J0+Ji+1,1):
            sum+=iter_count(J_list,M_values,i+1,J=Jf)
        return sum
    else:
        if M_values == 'pos':
            for Jf in np.arange(abs(J0-Ji),J0+Ji+1,1):
                sum+=Jf+1
            return sum
        elif M_values == 'none':
            for Jf in np.arange(abs(J0-Ji),J0+Ji+1,1):
                sum+=1
            return sum

def

#
#     if len(J_list)
#         num = 0
#         prevJ = J_list[0]
#         addJ = J_list[1:]
#         for i in range(len(addJ)):
#             newJ = addJ[i]
#             sumJ = np.arange(abs(prevJ-newJ),prevJ+newJ+1,1)
#
# def add_J(J1,J2):
#     if J1 is
#     size = len(np.arange(abs(J1-J2),J1+J2+1,1))
#
#
#
#      for Ji in ang_mom_list[1:]:
#
# def loop_rec(y, n):
#     if n >= 1:
#         for x in range(y):
#             loop_rec(y, n - 1)
#     else:
#        whatever()
#
# def loop_angmom(Jlist)
#
# def count_M(M_values):
#     if M_values==''
#
# else:
#     if

def q_numbers_bBJ_new(N_range,Lambda,S=1/2,I_list=[0,1/2],M_values='all'):
    IYb=I_list[0]
    iH = I_list[-1]
    Nmin,Nmax=N_range[0],N_range[-1]
    if Nmin<abs(Lambda):
        print('Nmin must be >= L')
        Nmin=abs(Lambda)
    if IYb == 0 or iH == 0:
        q_str = ['L','N','J','F','M']
        I = max(IYb,iH)
        q_numbers = {}
        for q in q_str:
            q_numbers[q] = []
        for N in np.arange(Nmin,Nmax+1,1):
            for J in np.arange(abs(N-S),abs(N+S)+1,1):
                for F in np.arange(abs(J-I),abs(J+I)+1,1):
                    if M_values=='none':
                        M=abs(F) % 1
                        for L in {True:[0], False:[-Lambda,Lambda]}[Lambda==0]:
                            values = [L,N,J,F,M]
                            for q,val in zip(q_str,values):
                                q_numbers[q].append(val+0)    #looks weird but adding 0 converts -0 to 0
                    else:
                        if M_values=='all':
                            Mmin = -F
                        elif M_values=='pos':
                            Mmin = abs(F) % 1
                        for M in np.arange(Mmin,F+1,1):
                            for L in {True:[0], False:[-Lambda,Lambda]}[Lambda==0]:
                                values = [L,N,J,F,M]
                                for q,val in zip(q_str,values):
                                    q_numbers[q].append(val+0)    #looks weird but adding 0 converts -0 to 0
    else:
        q_str = ['L','N','J','F1','F','M']
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
    return q_numbers

def q_numbers_bBJ_old(N_range,Lambda,S=1/2,I_list=[0,1/2],M_values='all'):
    IYb=I_list[0]
    iH = I_list[-1]
    Nmin,Nmax=N_range[0],N_range[-1]
    if Nmin<abs(Lambda):
        print('Nmin must be >= L')
        Nmin=abs(Lambda)
    if IYb == 0 or iH == 0:
        q_str = ['L','N','J','F','M']
        I = max(IYb,iH)
        q_numbers = {}
        for q in q_str:
            q_numbers[q] = []
        for N in np.arange(Nmin,Nmax+1,1):
            for J in np.arange(abs(N-S),abs(N+S)+1,1):
                for F in np.arange(abs(J-I),abs(J+I)+1,1):
                    if M_values=='none':
                        M=abs(F) % 1
                        for L in {True:[0], False:[-Lambda,Lambda]}[Lambda==0]:
                            values = [L,N,J,F,M]
                            for q,val in zip(q_str,values):
                                q_numbers[q].append(val+0)    #looks weird but adding 0 converts -0 to 0
                    else:
                        if M_values=='all':
                            Mmin = -F
                        elif M_values=='pos':
                            Mmin = abs(F) % 1
                        for M in np.arange(Mmin,F+1,1):
                            for L in {True:[0], False:[-Lambda,Lambda]}[Lambda==0]:
                                values = [L,N,J,F,M]
                                for q,val in zip(q_str,values):
                                    q_numbers[q].append(val+0)    #looks weird but adding 0 converts -0 to 0
    else:
        q_str = ['L','N','J','F1','F','M']
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
