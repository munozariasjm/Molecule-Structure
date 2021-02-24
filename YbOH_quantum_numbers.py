import numpy as np


#########################   DEVELOPMENT  BEGINNING   ###################################
def count_momenta(J_list, M_values):
    if M_values == 'all':
         num = 1
         for J in J_list:
             num*=(2*J+1)
    else:
        J_coupled = recursive_add_J(J_list)
        if M_values == 'pos':
            num = np.sum(J_coupled + 1 - (J_coupled[-1] % 1))
        elif M_values == 'none':
            num = len(J_coupled)
    return int(num)

def recursive_add_J(J_list):
    n = len(J_list)
    if n==2:
        J_sum = list_add_J(*J_list)
        return J_sum
    elif n>2:
        J_main = J_list[:-1]
        J_add = J_list[-1]
        return list_add_J(recursive_add_J(J_main),J_add)

def base_add_J(J1,J2):
    J_sum = np.arange(abs(J1-J2),J1+J2+1, 1)
    return J_sum

def list_add_J(J1,J2): #J2 cannot be a list
    if isinstance(J1, (list, tuple, np.ndarray)) and not isinstance(J2, (list, tuple, np.ndarray)):
        J_sum = []
        for _J in J1:
            J_sum.extend(base_add_J(_J,J2))
        return np.sort(np.array(J_sum))
    elif not isinstance(J1, (list, tuple, np.ndarray)) and not isinstance(J2, (list, tuple, np.ndarray)):
        return base_add_J(J1,J2)
    elif isinstance(J2, (list, tuple, np.ndarray)):
        if not isinstance(J1, (list, tuple, np.ndarray)):
            return list_add_J(J2,J1)
        else:
            print('J1 and J2 cannot both be lists')
            return None

def q_numbers_bBJ_new(N_range,Lambda,S=1/2,I_list=[0,1/2],M_values='all'):
    IYb=I_list[0]
    iH = I_list[-1]
    Nmin,Nmax=N_range[0],N_range[-1]
    if Nmin<abs(Lambda):
        print('Nmin must be >= L')
        Nmin=abs(Lambda)
    dim = 0
    for N in np.arange(Nmin,Nmax+1,1):
        dim+=count_momenta([N,S,*I_list],M_values)
    if Lambda != 0:
        dim*=2
    if IYb == 0 or iH == 0:
        q_str = ['L','N','J','F','M']
        I = max(IYb,iH)
        q_numbers = {}
        for q in q_str:
            q_numbers[q] = np.zeros(dim)
        i=0
        for N in np.arange(Nmin,Nmax+1,1):
            for J in np.arange(abs(N-S),abs(N+S)+1,1):
                for F in np.arange(abs(J-I),abs(J+I)+1,1):
                    if M_values=='none':
                        M=abs(F) % 1
                        for L in {True:[0], False:[-Lambda,Lambda]}[Lambda==0]:
                            values = [L,N,J,F,M]
                            for q,val in zip(q_str,values):
                                q_numbers[q][i] = val+0    #looks weird but adding 0 converts -0 to 0
                            i+=1
                    else:
                        if M_values=='all':
                            Mmin = -F
                        elif M_values=='pos':
                            Mmin = abs(F) % 1
                        for M in np.arange(Mmin,F+1,1):
                            for L in {True:[0], False:[-Lambda,Lambda]}[Lambda==0]:
                                values = [L,N,J,F,M]
                                for q,val in zip(q_str,values):
                                    q_numbers[q][i] = val+0    #looks weird but adding 0 converts -0 to 0
                                i+=1
    else:
        q_str = ['L','N','J','F1','F','M']
        q_numbers = {}
        for q in q_str:
            q_numbers[q] = np.zeros(dim)
        i=0
        for N in np.arange(Nmin,Nmax+1,1):
            for J in np.arange(abs(N-S),abs(N+S)+1,1):
                for F1 in np.arange(abs(J-IYb),abs(J+IYb)+1,1):
                    for F in np.arange(abs(F1-iH),abs(F1+iH)+1,1):
                        if M_values=='none':
                            M=abs(F) % 1
                            for L in {True:[0], False:[-Lambda,Lambda]}[Lambda==0]:
                                values = [L,N,J,F1,F,M]
                                for q,val in zip(q_str,values):
                                    q_numbers[q][i] = val+0    #looks weird but adding 0 converts -0 to 0
                                i+=1
                        else:
                            if M_values=='all':
                                Mmin = -F
                            elif M_values=='pos':
                                Mmin = abs(F) % 1
                            for M in np.arange(Mmin,F+1,1):
                                for L in {True:[0], False:[-Lambda,Lambda]}[Lambda==0]:
                                    values = [L,N,J,F1,F,M]
                                    for q,val in zip(q_str,values):
                                        q_numbers[q][i] = val+0    #looks weird but adding 0 converts -0 to 0
                                    i+=1
    return q_numbers


#########################   DEVELOPMENT  END   ###################################


def q_numbers_bBJ(N_range,Lambda,S=1/2,I_list=[0,1/2],M_values='all'):
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
        for N in np.arange(Nmin,Nmax+1,1):
            for J in np.arange(abs(N-S),abs(N+S)+1,1):
                for F1 in np.arange(abs(J-IYb),abs(J+IYb)+1,1):
                    for F in np.arange(abs(F1-iH),abs(F1+iH)+1,1):
                        if M_values=='none':
                            M=abs(F) % 1
                            for L in {True:[0], False:[-Lambda,Lambda]}[Lambda==0]:
                                values = [L,N,J,F1,F,M]
                                for q,val in zip(q_str,values):
                                    q_numbers[q].append(val+0)    #looks weird but adding 0 converts -0 to 0
                        else:
                            if M_values=='all':
                                Mmin = -F
                            elif M_values=='pos':
                                Mmin = abs(F) % 1
                            for M in np.arange(Mmin,F+1,1):
                                for L in {True:[0], False:[-Lambda,Lambda]}[Lambda==0]:
                                    values = [L,N,J,F1,F,M]
                                    for q,val in zip(q_str,values):
                                        q_numbers[q].append(val+0)    #looks weird but adding 0 converts -0 to 0
    return q_numbers


def q_numbers_aBJ(N_range,Lambda=1,S=1/2,I_list=[0,1/2],M_values='all',Omega_values=[1/2]):
    IYb=I_list[0]
    iH = I_list[-1]
    Nmin,Nmax=N_range[0],N_range[-1]
    Jmin = abs(Nmin-S)
    Jmax = abs(Nmax+S)
    q_numbers = {}
    for q in q_str:
        q_numbers[q] = []
    if Nmin<abs(Lambda):
        print('Nmin must be >= L')
        Nmin=abs(Lambda)
    if IYb == 0 or iH == 0:
        q_str = ['L','Sigma','Omega','J','F','M']
        I = max(IYb,iH)
        for J in np.arange(Jmin,Jmax+1,1):
            for F in np.arange(abs(J-I),abs(J+I)+1,1):
                if M_values=='none':
                    M=abs(F) % 1
                    for Sigma in np.arange(-abs(S),abs(S)+1,1):
                        for L in {True:[0], False:[-Lambda,Lambda]}[Lambda==0]:
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
                            for L in {True:[0], False:[-Lambda,Lambda]}[Lambda==0]:
                                Omega=L+Sigma
                                if abs(Omega) not in Omega_values:
                                    continue
                                else:
                                    values = [L,Sigma,Omega,J,F,M]
                                for q,val in zip(q_str,values):
                                    q_numbers[q].append(val+0)    #looks weird but adding 0 converts -0 to 0
    else:
        q_str = ['L','Sigma','Omega','J','F1','F','M']
        I = max(IYb,iH)
        for J in np.arange(Jmin,Jmax+1,1):
            for F1 in np.arange(abs(J-IYb),abs(J+IYb)+1,1):
                for F in np.arange(abs(F1-iH),abs(F1+iH)+1,1):
                    if M_values=='none':
                        M=abs(F) % 1
                        for Sigma in np.arange(-abs(S),abs(S)+1,1):
                            for L in {True:[0], False:[-Lambda,Lambda]}[Lambda==0]:
                                Omega=L+Sigma
                                if abs(Omega) not in Omega_values:
                                    continue
                                else:
                                    values = [L,Sigma,Omega,J,F1,F,M]
                                for q,val in zip(q_str,values):
                                    q_numbers[q].append(val+0)    #looks weird but adding 0 converts -0 to 0
                    else:
                        if M_values=='all':
                            Mmin = -F
                        elif M_values=='pos':
                            Mmin = abs(F) % 1
                        for M in np.arange(Mmin,F+1,1):
                            for Sigma in np.arange(-abs(S),abs(S)+1,1):
                                for L in {True:[0], False:[-Lambda,Lambda]}[Lambda==0]:
                                    Omega=L+Sigma
                                    if abs(Omega) not in Omega_values:
                                        continue
                                    else:
                                        values = [L,Sigma,Omega,J,F1,F,M]
                                    for q,val in zip(q_str,values):
                                        q_numbers[q].append(val+0)    #looks weird but adding 0 converts -0 to 0
    return q_numbers

def q_numbers_decoupled(N_range,Lambda=0, S=1/2, I_list=[0,1/2],M_values='all'):
    IYb=I_list[0]
    iH = I_list[-1]
    if M_values == 'none':
        print('Cannot construct decoupled basis without M values')
        return {}
    if Nmin<abs(Lambda):
        print('Nmin must be >= L')
        Nmin=abs(Lambda)
    elif IYb == 0 or iH == 0:
        Nmin,Nmax=N_range[0],N_range[-1]
        I = max(IYb,iH)
        q_str = ['L','N','M_N','M_S','M_I','M_F']
        q_numbers = {}
        for q in q_str:
            q_numbers[q] = []
        for L in {True:[0], False:[-Lambda,Lambda]}[Lambda==0]:
            for N in np.arange(Nmin,Nmax+1,1):
                for M_N in np.arange(-abs(N),abs(N)+1,1):
                    for M_S in np.arange(-abs(S),abs(S)+1,1):
                        for M_I in np.arange(-abs(I),abs(I)+1,1):
                            M_F = M_N + M_S + M_I
                            if (M_values == 'pos' and M_F<0):
                                continue
                            else:
                                values = [L,N,M_N,M_S,M_I,M_F]
                                for q,val in zip(q_str,values):
                                    q_numbers[q].append(val+0)    #looks weird but adding 0 converts -0 to 0
    else:
        for L in {True:[0], False:[-Lambda,Lambda]}[Lambda==0]:
            for N in np.arange(Nmin,Nmax+1,1):
                for M_N in np.arange(-abs(N),abs(N)+1,1):
                    for M_S in np.arange(-abs(S),abs(S)+1,1):
                        for M_IYb in np.arange(-abs(IYb),abs(IYb)+1,1):
                            for M_iH in np.arange(-abs(IYb),abs(IYb)+1,1):
                            M_F = M_N + M_S + M_IYb + M_iH
                            if (M_values=='pos' and M_F<0):
                                continue
                            else:
                                for L in range(-Lambda,Lambda+1,2):
                                    values = [L,N,M_N,M_S,M_IYb,M_iH,M_F]
                                    for q,val in zip(q_str,values):
                                        q_numbers[q].append(val+0)    #looks weird but adding 0 converts -0 to 0
    return q_numbers

def q_numbers_bBS(N_range,Lambda,S=1/2,I_list=[5/2,1/2],M_values='all'):
    IYb=I_list[0]
    iH = I_list[-1]
    Nmin,Nmax=N_range[0],N_range[-1]
    if Nmin<abs(Lambda):
        print('Nmin must be >= L')
        Nmin=abs(Lambda)
    if IYb == 0 or iH == 0:
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
