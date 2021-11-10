import sympy as sy
import numpy as np
from sympy.physics.wigner import wigner_3j,wigner_6j,wigner_9j

########## Matrix Elements for YbOH ##############


'''
Each matrix element is a function of an input state vector and an output
state vector.
The matrix elements are packaged into dictionaries. Each Hund's Case in the
molecule gets a separate dictionary. Dictionaries are at the end of the file.
'''



def kronecker(a,b):         # Kronecker delta function
    if a==b:
        return 1
    else:
        return 0

def b2a_matrix(a,b,S=1/2):
    if not kronecker(a['L'],b['L'])*kronecker(a['J'],b['J'])*kronecker(a['F'],b['F'])*kronecker(a['M'],b['M']):
        return 0
    else:
        if 'F1' in b.keys():
            if not kronecker(a['F1'],b['F1']):
                return 0
        return (-1)**(b['N']-S+a['Omega'])*np.sqrt(2*b['N']+1)*wigner_3j(a['J'],S,b['N'],a['Omega'],-a['Sigma'],-a['L'])

def decouple_b_174(dcpl,b,S=1/2,I=1/2): #dcpl = decoupled
    if not kronecker(dcpl['M_F'],b['M'])*kronecker(dcpl['L'],b['L'])*kronecker(dcpl['N'],b['N']):
        return 0
    else:
        M_J=dcpl['M_N']+dcpl['M_S']
        return (-1)**(I-b['J']+dcpl['M_F']+S-b['N']+M_J)*np.sqrt((2*b['F']+1)*(2*b['J']+1))*\
            wigner_3j(b['J'],I,b['F'],M_J,dcpl['M_I'],-dcpl['M_F'])*wigner_3j(b['N'],S,b['J'],dcpl['M_N'],dcpl['M_S'],-M_J)

def bBS_2_bBJ_matrix(bBS, bBJ, S=1/2, I = 5/2):
    if 'F1' in bBJ.keys():
        F = bBJ['F1']
        if not kronecker(bBS['F1'],bBJ['F1'])*kronecker(bBS['M'],bBJ['M'])*kronecker(bBS['N'],bBJ['N'])*kronecker(bBS['F'],bBJ['F']):
            return 0
    else:
        F = bBJ['F']
        if not kronecker(bBS['M'],bBJ['M'])*kronecker(bBS['N'],bBJ['N'])*kronecker(bBS['F'],bBJ['F']):
            return 0
    N = bBJ['N']
    G = bBS['G']
    J = bBJ['J']
    return (-1)**(I+S+F+N)*np.sqrt((2*G+1)*(2*J+1))*wigner_6j(I,S,G,N,F,J)


########## Case bBJ ##############

def Rot_bBJ(L0,N0,J0,F0,M0,L1,N1,J1,F1,M1,S=1/2,I=1/2):
    if not kronecker(L0,L1)*kronecker(N0,N1)*kronecker(F0,F1)*kronecker(M0,M1)*kronecker(J0,J1):
        return 0
    else:
        return N0*(N0+1)-L0**2

def SR_bBJ(L0,N0,J0,F0,M0,L1,N1,J1,F1,M1,S=1/2,I=1/2):
    if not kronecker(L0,L1)*kronecker(J0,J1)*kronecker(N0,N1)*kronecker(F0,F1)*kronecker(M0,M1):
        return 0
    else:
        return (-1)**(N0+J0+S)*np.sqrt(S*(S+1)*(2*S+1)*N0*(N0+1)*(2*N0+1))*\
            wigner_6j(S,N0,J0,N0,S,1)

def IS_bBJ(L0,N0,J0,F0,M0,L1,N1,J1,F1,M1,S=1/2,I=1/2):
    if not kronecker(L0,L1)*kronecker(N0,N1)*kronecker(F0,F1)*kronecker(M0,M1):
        return 0
    else:
        return (-1)**(J1 + F0 + I + J0 + N0 + S + 1)*\
            np.sqrt((2*J1+1)*(2*J0+1)*S*(S+1)*(2*S+1)*I*(I+1)*(2*I+1))*\
            wigner_6j(I,J1,F0,J0,I,1)*wigner_6j(J0,S,N0,S,J1,1)

def T2IS_bBJ(L0,N0,J0,F0,M0,L1,N1,J1,F1,M1,S=1/2,I=1/2):
    if not kronecker(F0,F1)*kronecker(M0,M1)*kronecker(L0,L1):
        return 0
    else:
        return -np.sqrt(5/3)*(-1)**(3*F0-2*M0+I+J0+N0-L0)*np.sqrt((2*I+1)*(I+1)*I)*\
            np.sqrt((2*J0+1)*(2*J1+1)*3*(2*S+1)*(S+1)*S)*wigner_6j(I,J1,F0,J0,I,1)*\
            wigner_9j(S,N1,J1,1,2,1,S,N0,J0)*np.sqrt((2*N0+1)*(2*N1+1))*wigner_3j(N0,2,N1,-L0,0,L1)


def Iz_bBJ(L0,N0,J0,F0,M0,L1,N1,J1,F1,M1,S=1/2,I=1/2):
    if not kronecker(F0,F1)*kronecker(M0,M1):
        return 0
    else:
        return (-1)**(I+F0+S-L0+1)*np.sqrt((2*J1+1)*(2*J0+1)*(2*N0+1)*(2*N1+1)*I*(I+1)*(2*I+1))*\
            wigner_6j(J1,I,F0,I,J0,1)*wigner_6j(N1,J1,S,J0,N0,1)*wigner_3j(N0,1,N1,-L0,0,L1)


def Sz_bBJ(L0,N0,J0,F0,M0,L1,N1,J1,F1,M1,S=1/2,I=1/2):
    if not kronecker(F0,F1)*kronecker(M0,M1)*kronecker(J0,J1):
        return 0
    else:
        return (-1)**(N1+N0+S+J0-L0)*np.sqrt((2*N0+1)*(2*N1+1)*S*(S+1)*(2*S+1))*\
            wigner_6j(N1,S,J0,S,N0,1)*wigner_3j(N0,1,N1,-L0,0,L1)

def NzSz_bBJ(L0,N0,J0,F0,M0,L1,N1,J1,F1,M1,S=1/2,I=1/2):
    if not kronecker(F0,F1)*kronecker(M0,M1)*kronecker(J0,J1)*kronecker(L0,L1):
        return 0
    else:
        return L0*(-1)**(N1+N0+S+J0-L0)*np.sqrt((2*N0+1)*(2*N1+1)*S*(S+1)*(2*S+1))*\
            wigner_6j(N1,S,J0,S,N0,1)*wigner_3j(N0,1,N1,-L0,0,L1)

def ZeemanZ_bBJ(L0,N0,J0,F0,M0,L1,N1,J1,F1,M1,S=1/2,I=1/2):
    if not kronecker(L0,L1)*kronecker(M0,M1)*kronecker(N0,N1):
        return 0
    else:
        return (-1)**(F0-M0+F1+2*J0+I+N0+S)*np.sqrt((2*F0+1)*(2*F1+1)*(2*J0+1)*(2*J1+1)*S*(S+1)*(2*S+1))*\
            wigner_6j(F0,J0,I,J1,F1,1)*wigner_6j(J0,S,N0,S,J1,1)*wigner_3j(F0,1,F1,-M0,0,M1)

#def ZeemanPlus_bBJ():

#def ZeemanMinus_bBJ():


def StarkZ_bBJ(L0,N0,J0,F0,M0,L1,N1,J1,F1,M1,S=1/2,I=1/2):
    if not kronecker(M0,M1):
        return 0
    else:
        return (-1)**(F0-M0+J0+I+2+F1+S+J1-L0)*np.sqrt((2*F0+1)*(2*F1+1)*(2*J0+1)*(2*J1+1)*(2*N0+1)*(2*N1+1))*\
            wigner_3j(F0,1,F1,-M0,0,M1)*wigner_6j(J1,F1,I,F0,J0,1)*wigner_6j(N1,J1,S,J0,N0,1)*wigner_3j(N0,1,N1,-L0,0,L1)

def lD_bBJ(L0,N0,J0,F0,M0,L1,N1,J1,F1,M1,S=1/2,I=1/2):
    if not kronecker(N0,N1)*kronecker(M0,M1)*kronecker(J0,J1)*kronecker(F0,F1)*(not kronecker(L0,L1)):
        return 0
    else:
        delta_L = L0-L1
        return N0*(N0+1)*(kronecker(delta_L,-2)*(-1)**(-abs(L0)) + kronecker(delta_L,2)*(-1)**(abs(L0)))


########## Case bBS ##############

def Rot_bBS(L0,N0,G0,F10,F0,M0,L1,N1,G1,F11,F1,M1,S=1/2,I=5/2,iH=1/2):
    if not kronecker(L0,L1)*kronecker(F0,F1)*kronecker(F10, F11)*kronecker(M0,M1)*kronecker(G0,G1)*kronecker(N0,N1):
        return 0
    else:
        return (-1)**(-2*M0+2*iH+2*2*G0)*N0*(N0+1)-L0**2

def SR_bBS(L0,N0,G0,F10,F0,M0,L1,N1,G1,F11,F1,M1,S=1/2,I=5/2,iH=1/2):
    if not kronecker(L0,L1)*kronecker(F0,F1)*kronecker(F10, F11)*kronecker(M0,M1)*kronecker(N0,N1):
        return 0
    else:
        return (-1)**(-2*M0+2*iH+3*F10+N1+G0+G1+I+S+1)*np.sqrt((2*G0+1)*(2*G1+1)*N0*(N0+1)*(2*N0+1)*S*(S+1)*(2*S+1))*\
            wigner_6j(N1,G1,F10,G0,N0,1)*wigner_6j(S,G1,I,G0,S,1)

def NzSz_bBS(L0,N0,G0,F10,F0,M0,L1,N1,G1,F11,F1,M1,S=1/2,I=5/2,iH=1/2):
    if not kronecker(L0,L1)*kronecker(F0,F1)*kronecker(F10, F11)*kronecker(M0,M1)*kronecker(N0,N1):
        return 0
    else:
        return L0*(-1)**(-2*M0+2*iH+2*F10+N1+G0+F10+N0-L0+G1+S+I+1)*np.sqrt((2*N0+1)*(2*N1+1)*(2*G0+1)*(2*G1+1))*\
            wigner_6j(N1,G1,F10,G0,N0,1)*wigner_6j(S,G1,I,G0,S,1)*wigner_3j(N0,1,N1,-L0,0,L1)*np.sqrt(S*(S+1)*(2*S+1))

def ISYb_bBS(L0,N0,G0,F10,F0,M0,L1,N1,G1,F11,F1,M1,S=1/2,I=5/2,iH=1/2):
    if not kronecker(L0,L1)*kronecker(F0,F1)*kronecker(F10, F11)*kronecker(M0,M1)*kronecker(N0,N1)*kronecker(G0,G1):
        return 0
    else:
        return (-1)**(-2*M0+2*iH+2*N0+2*G0+I+S+G0)*np.sqrt(S*(S+1)*(2*S+1)*I*(I+1)*(2*I+1))*\
            wigner_6j(I,S,G0,S,I,1)

def T2ISYb_bBS(L0,N0,G0,F10,F0,M0,L1,N1,G1,F11,F1,M1,S=1/2,I=5/2,iH=1/2):
    if not kronecker(F0,F1)*kronecker(F10, F11)*kronecker(M0,M1)*kronecker(L0,L1):
        return 0
    else:
        return (-1)**(-2*M0+2*iH+3*F10+N0+G1+N0-L0)*wigner_9j(G0,G1,2,I,I,1,S,S,1)*\
            np.sqrt((2*N0+1)*(2*N1+1)*5*(2*G0+1)*(2*G1+1)*I*(I+1)*(2*I+1)*S*(S+1)*(2*S+1))*\
            wigner_6j(N1,G1,F10,G0,N0,2)*wigner_3j(N0,2,N1,-L0,0,L1)

def T2QYb_bBS(L0,N0,G0,F10,F0,M0,L1,N1,G1,F11,F1,M1,S=1/2,I=5/2,iH=1/2):
    if not kronecker(F0,F1)*kronecker(F10, F11)*kronecker(M0,M1)*kronecker(L0,L1):
        return 0
    else:
        return (-1)**(-2*M0+2*iH+3*F10+N0+G1+N0-L0+G1+I+S+2)*np.sqrt((2*N0+1)*(2*N1+1)*(2*G0+1)*(2*G1+1))*\
            I*(2*I-1)/np.sqrt(6)/wigner_3j(I,2,I,-I,0,I)*wigner_3j(N0,2,N1,-L0,0,L1)*\
            wigner_6j(N1,G1,F10,G0,N0,2)*wigner_6j(I,G1,S,G0,I,2)

def ISH_bBS(L0,N0,G0,F10,F0,M0,L1,N1,G1,F11,F1,M1,S=1/2,I=5/2,iH=1/2):
    if not kronecker(F0,F1)*kronecker(N0, N1)*kronecker(M0,M1)*kronecker(L0,L1):
        return 0
    else:
        return (-1)**(2*F0-2*M0+iH+F0+F10+F11+N0+G0+1+G1+I+S+1)*\
        np.sqrt(iH*(iH+1)*(2*iH+1)*(2*F10+1)*(2*F11+1)*(2*G0+1)*(2*G1+1)*S*(S+1)*(2*S+1))*\
        wigner_6j(iH,F11,F0,F10,iH,1)*wigner_6j(G1,F11,N0,F10,G0,1)*wigner_6j(S,G1,I,G0,S,1)

def T2ISH_bBS(L0,N0,G0,F10,F0,M0,L1,N1,G1,F11,F1,M1,S=1/2,I=5/2,iH=1/2):
    if not kronecker(F0,F1)*kronecker(M0,M1)*kronecker(L0,L1):
        return 0
    else:
        return (-1)**(3*F0-2*M0+iH+F10+N0-L0+G1+S+I+1)*\
            np.sqrt(iH*(iH+1)*(2*iH+1)*(2*F10+1)*(2*F11+1)*3*(2*G0+1)*(2*G1+1)*S*(S+1)*(2*S+1))*\
            wigner_6j(iH,F11,F0,F10,iH,1)*wigner_9j(F10,F11,1,N0,N1,2,G0,G1,1)*\
            wigner_3j(N0,2,N1,-L0,0,L1)*wigner_6j(S,G1,I,G0,S,1)

def StarkZ_bBS(L0,N0,G0,F10,F0,M0,L1,N1,G1,F11,F1,M1,S=1/2,I=5/2,iH=1/2):
    if not kronecker(G0,G1)*kronecker(M0,M1)*kronecker(L0,L1):
        return 0
    else:
        return (-1)**(F0-M0+F1+F10+iH+F11+2*N0+G0+2-L0)*wigner_3j(F0,1,F1,-M0,0,M1)*\
            np.sqrt((2*F0+1)*(2*F1+1)*(2*F10+1)*(2*F11+1)*(2*N0+1)*(2*N1+1))*\
            wigner_6j(F11,F1,iH,F0,F10,1)*wigner_6j(N1,F11,G0,F10,N0,1)*wigner_3j(N0,1,N1,-L0,0,L1)

def ZeemanZ_bBS(L0,N0,G0,F10,F0,M0,L1,N1,G1,F11,F1,M1,S=1/2,I=5/2,iH=1/2):
    if not kronecker(N0,N1)*kronecker(M0,M1)*kronecker(L0,L1):
        return 0
    else:
        return (-1)**(F0-M0+F1+F10+iH+1+F11+N0+G0+1+G0+I+S+1)*wigner_3j(F0,1,F1,-M0,0,M1)*\
            np.sqrt((2*F0+1)*(2*F1+1)*(2*F10+1)*(2*F11+1)*(2*G0+1)*(2*G1+1)*S*(S+1)*(2*S+1))*\
            wigner_6j(F11,F1,iH,F0,F10,1)*wigner_6j(G1,F11,N0,F10,G0,1)*wigner_6j(S,G1,I,G0,S,1)

def lD_bBS(L0,N0,G0,F10,F0,M0,L1,N1,G1,F11,F1,M1,S=1/2,I=5/2,iH=1/2):
    if not kronecker(N0,N1)*kronecker(M0,M1)*kronecker(G0,G1)*kronecker(F10,F11)*kronecker(F0,F1)*(not kronecker(L0,L1)):
        return 0
    else:
        delta_L = L0-L1
        return N0*(N0+1)*(kronecker(delta_L,-2)*(-1)**(-abs(L0)) + kronecker(delta_L,2)*(-1)**(abs(L0)))


#Using combined tensor form

def MQM_bBS_old(L0,N0,G0,F10,F0,M0,L1,N1,G1,F11,F1,M1,S=1/2,I=5/2,iH=1/2):
    if not kronecker(L0,L1)*kronecker(F0,F1)*kronecker(F10,F11)*kronecker(M0,M1):
        return 0
    else:
        return (-1)**(-2*M0+2*iH+2*F10+N0-L0)*np.sqrt((2*G0+1)*(2*G1+1)*(3)*(2*F11+1))*\
            wigner_9j(F10,F11,0,G0,G1,1,N0,N1,1)*wigner_9j(G0,G1,1,I,I,2,S,S,1)*\
            np.sqrt((2*N0+1)*(2*N1+1)*S*(S+1)*(2*S+1))*I*(2*I-1)/(np.sqrt(6)*wigner_3j(I,2,I,-I,0,I))*\
            wigner_3j(N0,1,N1,-L0,0,L1)

def MQM_bBS(L0,N0,G0,F10,F0,M0,L1,N1,G1,F11,F1,M1,S=1/2,I=5/2,iH=1/2):
    if not kronecker(L0,L1)*kronecker(F0,F1)*kronecker(F10,F11)*kronecker(M0,M1):
        return 0
    else:
        return (-1)**(-2*M0+2*iH+2*F10+N1+G0+F10)*wigner_6j(N1,G1,F10,G0,N0,1)*\
            (-1)**(N0-L0)*wigner_3j(N0,1,N1,-L0,0,L1)*np.sqrt((2*N0+1)*(2*N1+1))*\
            np.sqrt((2*G0+1)*(2*G1+1)*3)*wigner_9j(G0,G1,1,I,I,2,S,S,1)*\
            np.sqrt((2*S+1)*(S+1)*S)*I*(2*I-1)/np.sqrt(6)/wigner_3j(I,2,I,-I,0,I)
#Using Sz(3Iz-I2) format



def EDM_bBS(L0,N0,G0,F10,F0,M0,L1,N1,G1,F11,F1,M1,S=1/2,I=5/2,iH=1/2):
    if not kronecker(L0,L1)*kronecker(F0,F1)*kronecker(F10,F11)*kronecker(M0,M1):
        return 0
    else:
        return (-1)**(-2*M0+2*iH+2*F10+N1+G0+F10+N0-L0+G1+S+I+1)*np.sqrt((2*N0+1)*(2*N1+1)*(2*G0+1)*(2*G1+1))*\
            wigner_6j(N1,G1,F10,G0,N0,1)*wigner_6j(S,G1,I,G0,S,1)*wigner_3j(N0,1,N1,-L0,0,L1)*np.sqrt(S*(S+1)*(2*S+1))







########## 174YbOH Case aBJ ##############

def Rot_174_aBJ(L0,Sigma0,Omega0,J0,F0,M0,L1,Sigma1,Omega1,J1,F1,M1,S=1/2,I=1/2):
    if not kronecker(L0,L1)*kronecker(F0,F1)*kronecker(M0,M1)*kronecker(J0,J1):
        return 0
    else:
        return kronecker(Sigma0,Sigma1)*kronecker(Omega0,Omega1)*(J0*(J0+1)+S*(S+1)-2*Omega0*Sigma0)-\
            2*(-1)**(J0-(Omega0)+S-Sigma0)*np.sqrt((2*J0+1)*J0*(J0+1)*(2*S+1)*S*(S+1))*\
            sum([wigner_3j(J0,1,J1,-Omega0,q,Omega1)*wigner_3j(S,1,S,-Sigma0,q,Sigma1) for q in [-1,1]])

            #Note: in N^2 formulation, you get -2*Omega*Sigma. Brown uses R^2 form which results in
            #a term -Omega^2 - Sigma^2. You can show for us that Omega^2 + Sigma^2 = Lambda^2 + 2 Omega*Sigma

def SO_174_aBJ(L0,Sigma0,Omega0,J0,F0,M0,L1,Sigma1,Omega1,J1,F1,M1,S=1/2,I=1/2):
    if not kronecker(L0,L1)*kronecker(F0,F1)*kronecker(M0,M1)*kronecker(J0,J1)*kronecker(Sigma0,Sigma1)*kronecker(Omega0,Omega1):
        return 0
    else:
        return L0*Sigma0

def IL_174_aBJ(L0,Sigma0,Omega0,J0,F0,M0,L1,Sigma1,Omega1,J1,F1,M1,S=1/2,I=1/2):
    if not kronecker(L0,L1)*kronecker(F0,F1)*kronecker(M0,M1):
        return 0
    else:
        return L0*(-1)**(2*F0-2*M0-J1+I+F0+J0-Omega0)*wigner_6j(J1,I,F0,I,J0,1)*\
            wigner_3j(J0,1,J1,-Omega0,0,Omega1)*np.sqrt((2*J0+1)*(2*J1+1)*(2*I+1)*(I+1)*I)

def IS_174_aBJ(L0,Sigma0,Omega0,J0,F0,M0,L1,Sigma1,Omega1,J1,F1,M1,S=1/2,I=1/2):
    if not kronecker(L0,L1)*kronecker(F0,F1)*kronecker(M0,M1):
        return 0
    else:
        return (-1)**(2*F0-2*M0+J1+I+F0+S-Sigma0+J0-Omega0)*wigner_6j(J1,I,F0,I,J0,1)*\
            np.sqrt((2*S+1)*(S+1)*S*(2*J0+1)*(2*J1+1)*I*(I+1)*(2*I+1))*\
            sum([wigner_3j(J0,1,J1,-Omega0,q,Omega1)*wigner_3j(S,1,S,-Sigma0,q,Sigma1) for q in [-1,0,1]])

def IzSz_174_aBJ(L0,Sigma0,Omega0,J0,F0,M0,L1,Sigma1,Omega1,J1,F1,M1,S=1/2,I=1/2):
    if not kronecker(L0,L1)*kronecker(F0,F1)*kronecker(M0,M1)*kronecker(Sigma0,Sigma1):
        return 0
    else:
        return Sigma0*(-1)**(2*F0-2*M0+J1+I+F0+J0-Omega0)*wigner_6j(J1,I,F0,I,J0,1)*\
            np.sqrt((2*J0+1)*(2*J1+1)*I*(I+1)*(2*I+1))*wigner_3j(J0,1,J1,-Omega0,0,Omega1)

def T2IS_174_aBJ(L0,Sigma0,Omega0,J0,F0,M0,L1,Sigma1,Omega1,J1,F1,M1,S=1/2,I=1/2):
    if not kronecker(L0,L1)*kronecker(F0,F1)*kronecker(M0,M1):
        return 0
    else:
        return (-1)**(2*F0-2*M0+J1+I+F0+J0-Omega0+S-Sigma0)*wigner_6j(J1,I,F0,I,J0,1)*\
            np.sqrt(5*(2*J0+1)*(2*J1+1)*S*(S+1)*(2*S+1)*I*(I+1)*(2*I+1))*\
            sum([(-1)**(q)*wigner_3j(1,1,2,q,-q,0)*wigner_3j(S,1,S,-Sigma0,-q,Sigma1)*\
            wigner_3j(J0,1,J1,-Omega0,-q,Omega1) for q in [-1,0,1]])

def ZeemanLZ_174_aBJ(L0,Sigma0,Omega0,J0,F0,M0,L1,Sigma1,Omega1,J1,F1,M1,S=1/2,I=1/2):
    if not kronecker(L0,L1)*kronecker(Sigma0,Sigma1)*kronecker(M0,M1)*kronecker(Omega0,Omega1):
        return 0
    else:
        return L0*(-1)**(F0-M0)*wigner_3j(F0,1,F1,-M0,0,M1)*\
            (-1)**(F1+J0+I+1)*np.sqrt((2*F0+1)*(2*F1+1))*wigner_6j(J1,F1,I,F0,J0,1)*\
            (-1)**(J0-Omega0)*np.sqrt((2*J0+1)*(2*J1+1))*wigner_3j(J0,1,J1,-Omega0,0,Omega1)

def ZeemanSZ_174_aBJ(L0,Sigma0,Omega0,J0,F0,M0,L1,Sigma1,Omega1,J1,F1,M1,S=1/2,I=1/2):
    if not kronecker(L0,L1)*kronecker(M0,M1):
        return 0
    else:
        return (-1)**(F0-M0)*wigner_3j(F0,1,F1,-M0,0,M1)*\
            (-1)**(F1+J0+I+1)*np.sqrt((2*F0+1)*(2*F1+1))*wigner_6j(J1,F1,I,F0,J0,1)*\
            (-1)**(J0-Omega0+S-Sigma0)*np.sqrt((2*J0+1)*(2*J1+1)*S*(S+1)*(2*S+1))*\
            sum([wigner_3j(J0,1,J1,-Omega0,q,Omega1)*wigner_3j(S,1,S,-Sigma0,q,Sigma1) for q in [-1,0,1]])

def ZeemanParityZ_174_aBJ(L0,Sigma0,Omega0,J0,F0,M0,L1,Sigma1,Omega1,J1,F1,M1,S=1/2,I=1/2):
    if kronecker(L0,L1)*(not kronecker(M0,M1)):
        return 0
    else:
        return (-1)*(-1)**(F0-M0)*wigner_3j(F0,1,F1,-M0,0,M1)*\
            (-1)**(F1+J0+I+1)*np.sqrt((2*F0+1)*(2*F1+1))*wigner_6j(J1,F1,I,F0,J0,1)*\
            np.sqrt((2*J0+1)*(2*J1+1)*S*(S+1)*(2*S+1))*\
            sum([kronecker(L1,L0-2*q)*(-1)**(J0-Omega0+S-Sigma0)*wigner_3j(J0,1,J1,-Omega0,q,Omega1)*wigner_3j(S,1,S,-Sigma0,-q,Sigma1) for q in [-1,1]])

def StarkZ_174_aBJ(L0,Sigma0,Omega0,J0,F0,M0,L1,Sigma1,Omega1,J1,F1,M1,S=1/2,I=1/2):
    if not kronecker(L0,L1)*kronecker(Sigma0,Sigma1):
        return 0
    else:
        return (-1)**(F0-M0)*wigner_3j(F0,1,F1,-M0,0,M1)*\
            (-1)**(F1+J0+I+1)*np.sqrt((2*F0+1)*(2*F1+1))*wigner_6j(J1,F1,I,F0,J0,1)*\
            (-1)**(J0-Omega0)*wigner_3j(J0,1,J1,-Omega0,0,Omega1)*np.sqrt((2*J0+1)*(2*J1+1))

def LambdaDoubling_174_aBJ(L0,Sigma0,Omega0,J0,F0,M0,L1,Sigma1,Omega1,J1,F1,M1,S=1/2,I=1/2):
    if not kronecker(F0,F1)*kronecker(M0,M1)*kronecker(J0,J1)*(not kronecker(L0,L1)):
        return 0
    else:
        return (-1)**(J0-Omega0+S-Sigma0+1)*np.sqrt((2*J0+1)*(J0+1)*J0*(2*S+1)*(S+1)*S)*\
            sum([kronecker(L1,L0+2*q)*wigner_3j(J0,1,J1,-Omega0,-q,Omega1)*wigner_3j(S,1,S,-Sigma0,q,Sigma1) for q in [-1,1]])

def TransitionDipole_174_aBJ(L0,Sigma0,Omega0,J0,F0,M0,L1,Sigma1,Omega1,J1,F1,M1,S=1/2,I=1/2):
    if not kronecker(Sigma0,Sigma1):
        return 0
    else:
        TDM_total = sum([(-1)**(p)*sum([(-1)**(F0-M0)*wigner_3j(F0,1,F1,-M0,p,M1)*(-1)**(F1+J0+I+1)*np.sqrt((2*F0+1)*(2*F1+1))*wigner_6j(J1,F1,I,F0,J0,1)*(-1)**(J0-Omega0)*np.sqrt((2*J0+1)*(2*J1+1))*wigner_3j(J0,1,J1,-Omega0,q,Omega1) for q in [-1,1]]) for p in range(-1,2)])
             # Since L is changing, can only get q=+-1 transitions
        # TDM_plus = sum([(-1)**(F0-M0)*wigner_3j(F0,1,F1,-M0,1,M1)*(-1)**(F1+J+I+1)*np.sqrt((2*F0+1)*(2*F1+1))*\
        #     wigner_6j(J1,F1,I,F0,J0,1)*(-1)**(J0-Omega0)*np.sqrt((2*J0+1)*(2*J1+1))*wigner_3j(J0,1,J1,-Omega0,q,Omega1) for q in range(-1,2)])
        # TDM_zero = sum([(-1)**(F0-M0)*wigner_3j(F0,1,F1,-M0,0,M1)*(-1)**(F1+J+I+1)*np.sqrt((2*F0+1)*(2*F1+1))*\
        #     wigner_6j(J1,F1,I,F0,J0,1)*(-1)**(J0-Omega0)*np.sqrt((2*J0+1)*(2*J1+1))*wigner_3j(J0,1,J1,-Omega0,q,Omega1) for q in range(-1,2)])
        # TDM_minus = sum([(-1)**(F0-M0)*wigner_3j(F0,1,F1,-M0,-1,M1)*(-1)**(F1+J+I+1)*np.sqrt((2*F0+1)*(2*F1+1))*\
        #     wigner_6j(J1,F1,I,F0,J0,1)*(-1)**(J0-Omega0)*np.sqrt((2*J0+1)*(2*J1+1))*wigner_3j(J0,1,J1,-Omega0,q,Omega1) for q in range(-1,2)])
        # TDM_total = TDM_plus + TDM_zero + TDM_minus
        return TDM_total

def TransitionDipole_174_aBJ_noM(L0,Sigma0,Omega0,J0,F0,M0,L1,Sigma1,Omega1,J1,F1,M1,S=1/2,I=1/2):
    if not kronecker(Sigma0,Sigma1):
        return 0
    elif not (kronecker(F0,F1) or kronecker(F0+1,F1) or kronecker(F0-1,F1)):
        return 0
    else:
        TDM_total = 1/np.sqrt((2*F1+1))*sum([
                        (-1)**(F1+J0+I+1)*np.sqrt((2*F0+1)*(2*F1+1))*\
                            wigner_6j(J1,F1,I,F0,J0,1)*(-1)**(J0-Omega0)*np.sqrt((2*J0+1)*(2*J1+1))*wigner_3j(J0,1,J1,-Omega0,q,Omega1)\
                            for q in [-1,1]]) # Since L is changing, can only get q=+-1 transitions
        # TDM_plus = sum([(-1)**(F0-M0)*wigner_3j(F0,1,F1,-M0,1,M1)*(-1)**(F1+J+I+1)*np.sqrt((2*F0+1)*(2*F1+1))*\
        #     wigner_6j(J1,F1,I,F0,J0,1)*(-1)**(J0-Omega0)*np.sqrt((2*J0+1)*(2*J1+1))*wigner_3j(J0,1,J1,-Omega0,q,Omega1) for q in range(-1,2)])
        # TDM_zero = sum([(-1)**(F0-M0)*wigner_3j(F0,1,F1,-M0,0,M1)*(-1)**(F1+J+I+1)*np.sqrt((2*F0+1)*(2*F1+1))*\
        #     wigner_6j(J1,F1,I,F0,J0,1)*(-1)**(J0-Omega0)*np.sqrt((2*J0+1)*(2*J1+1))*wigner_3j(J0,1,J1,-Omega0,q,Omega1) for q in range(-1,2)])
        # TDM_minus = sum([(-1)**(F0-M0)*wigner_3j(F0,1,F1,-M0,-1,M1)*(-1)**(F1+J+I+1)*np.sqrt((2*F0+1)*(2*F1+1))*\
        #     wigner_6j(J1,F1,I,F0,J0,1)*(-1)**(J0-Omega0)*np.sqrt((2*J0+1)*(2*J1+1))*wigner_3j(J0,1,J1,-Omega0,q,Omega1) for q in range(-1,2)])
        # TDM_total = TDM_plus + TDM_zero + TDM_minus
        return TDM_total

########## 173YbOH Case aBJ ##############

def SO_173_aBJ(L0,Sigma0,Omega0,J0,F10,F0,M0,L1,Sigma1,Omega1,J1,F11,F1,M1,S=1/2,I=5/2,iH=1/2):
    if not kronecker(L0,L1)*kronecker(F0,F1)*kronecker(M0,M1)*kronecker(J0,J1)*kronecker(Sigma0,Sigma1)*kronecker(Omega0,Omega1)*kronecker(F10,F11):
        return 0
    else:
        return L0*Sigma0

def Rot_173_aBJ(L0,Sigma0,Omega0,J0,F10,F0,M0,L1,Sigma1,Omega1,J1,F11,F1,M1,S=1/2,I=5/2,iH=1/2):
    if not kronecker(L0,L1)*kronecker(F0,F1)*kronecker(F10,F11)*kronecker(M0,M1)*kronecker(J0,J1)*kronecker(Sigma0,Sigma1)*kronecker(Omega0,Omega1):
        return 0
    else:
        return (J0*(J0+1)+S*(S+1)-Omega0**2 -Sigma0**2)-\
            2*(-1)**(J0-(Omega0)+S-Sigma0)*np.sqrt((2*J0+1)*J0*(J0+1)*(2*S+1)*S*(S+1))*\
            sum([wigner_3j(J0,1,J1,-Omega0,q,Omega1)*wigner_3j(S,1,S,-Sigma0,q,Sigma1) for q in [-1,1]])

def ILYb_173_aBJ(L0,Sigma0,Omega0,J0,F10,F0,M0,L1,Sigma1,Omega1,J1,F11,F1,M1,S=1/2,I=5/2,iH=1/2):
    if not kronecker(L0,L1)*kronecker(F0,F1)*kronecker(M0,M1)*kronecker(F10,F11)*kronecker(Sigma0,Sigma1)*kronecker(Omega0,Omega1):
        return 0
    else:
        return L0*(-1)**(J1+I+F10)*wigner_6j(J1,I,F10,I,J0,1)*\
            (-1)**(J0-Omega0)*wigner_3j(J0,1,J1,-Omega0,0,Omega1)*\
            np.sqrt((2*J0+1)*(2*J1+1)*(2*I+1)*(I+1)*I)

#Check derivation
def T2q2_ISYb_173_aBJ(L0,Sigma0,Omega0,J0,F10,F0,M0,L1,Sigma1,Omega1,J1,F11,F1,M1,S=1/2,I=5/2,iH=1/2):
    if not kronecker(F0,F1)*kronecker(M0,M1)*kronecker(F10,F11):
        return 0
    else:
        return (-1)**(J1+I+F10+J0-Omega0+S-Sigma0)*\
            wigner_6j(J1,I,F10,I,J0,1)*np.sqrt((2*J0+1)*(2*J1+1)*S*(S+1)*(2*S+1)*I*(I+1)*(2*I+1))*\
            sum([(-1)**(q)*kronecker(L0,-L1)*wigner_3j(S,1,S,-Sigma0,q,Sigma1)*wigner_3j(J0,1,J1,-Omega0,-q,Omega1) for q in [-1,1]])

#Check derivation
def T2q0_IYb_173_aBJ(L0,Sigma0,Omega0,J0,F10,F0,M0,L1,Sigma1,Omega1,J1,F11,F1,M1,S=1/2,I=5/2,iH=1/2):
    if not kronecker(F0,F1)*kronecker(M0,M1)*kronecker(F10,F11)*kronecker(L0,L1)*kronecker(Omega0,Omega1)*kronecker(L0,L1)*kronecker(Sigma0,Sigma1):
        return 0
    else:
        return 1/4*(-1)**(J1+I+F10+J0-Omega0)*wigner_6j(J1,I,F10,I,J0,2)*\
            wigner_3j(J0,2,J1,-Omega0,0,Omega1)*np.sqrt((2*J0+1)*(2*J1+1))/\
            wigner_3j(I,2,I,-I,0,I)

#Check derivation
def LambdaDoubling_173_aBJ(L0,Sigma0,Omega0,J0,F10,F0,M0,L1,Sigma1,Omega1,J1,F11,F1,M1,S=1/2,I=5/2,iH=1/2):
    if not kronecker(F0,F1)*kronecker(M0,M1)*kronecker(F10,F11)*kronecker(J0,J1)*(not kronecker(L0,L1)):
        return 0
    else:
        return (-1)**(J0-Omega0+S-Sigma0)*np.sqrt((2*J0+1)*(J0+1)*J0*(2*S+1)*(S+1)*S)*\
            sum([kronecker(L0,L1-2*q)*wigner_3j(J0,1,J1,-Omega0,-q,Omega1)*wigner_3j(S,1,S,-Sigma0,q,Sigma1) for q in [-1,1]])


def ZeemanLZ_173_aBJ(L0,Sigma0,Omega0,J0,F10,F0,M0,L1,Sigma1,Omega1,J1,F11,F1,M1,S=1/2,I=5/2,iH=1/2):
    if not kronecker(L0,L1)*kronecker(Sigma0,Sigma1):
        return 0
    else:
        return L0*(-1)**(F0-M0)*wigner_3j(F0,1,F1,-M0,0,M1)*\
            (-1)**(F1+F10+iH+1)*np.sqrt((2*F0+1)*(2*F1+1))*wigner_6j(F11,F1,iH,F0,F10,1)*\
            (-1)**(F11+J0+I+1)*np.sqrt((2*F10+1)*(2*F11+1))*wigner_6j(J1,F11,I,F10,J0,1)*\
            (-1)**(J0-Omega0)*np.sqrt((2*J0+1)*(2*J1+1))*wigner_3j(J0,1,J1,-Omega0,0,Omega1)

def ZeemanSZ_173_aBJ(L0,Sigma0,Omega0,J0,F10,F0,M0,L1,Sigma1,Omega1,J1,F11,F1,M1,S=1/2,I=5/2,iH=1/2):
    if not kronecker(L0,L1):
        return 0
    else:
        return (-1)**(F0-M0)*wigner_3j(F0,1,F1,-M0,0,M1)*\
            (-1)**(F1+F10+iH+1)*np.sqrt((2*F0+1)*(2*F1+1))*wigner_6j(F11,F1,iH,F0,F10,1)*\
            (-1)**(F11+J0+I+1)*np.sqrt((2*F10+1)*(2*F11+1))*wigner_6j(J1,F11,I,F10,J0,1)*\
            (-1)**(J0-Omega0+S-Sigma0)*np.sqrt((2*J0+1)*(2*J1+1)*S*(S+1)*(2*S+1))*\
            sum([wigner_3j(J0,1,J1,-Omega0,q,Omega1)*wigner_3j(S,1,S,-Sigma0,q,Sigma1) for q in [-1,0,1]])

def ZeemanParityZ_173_aBJ(L0,Sigma0,Omega0,J0,F10,F0,M0,L1,Sigma1,Omega1,J1,F11,F1,M1,S=1/2,I=5/2,iH=1/2):
    if kronecker(L0,L1):
        return 0
    else:
        return (-1)**(F0-M0)*wigner_3j(F0,1,F1,-M0,0,M1)*\
            (-1)**(F1+F10+iH+1)*np.sqrt((2*F0+1)*(2*F1+1))*wigner_6j(F11,F1,iH,F0,F10,1)*\
            (-1)**(F11+J0+I+1)*np.sqrt((2*F10+1)*(2*F11+1))*wigner_6j(J1,F11,I,F10,J0,1)*\
            np.sqrt((2*J0+1)*(2*J1+1)*S*(S+1)*(2*S+1))*\
            sum([kronecker(L0,L1-2*q)*(-1)**(J0-Omega0+S-Sigma0)*wigner_3j(J0,1,J1,-Omega0,-q,Omega1)*wigner_3j(S,1,S,-Sigma0,q,Sigma1) for q in [-1,1]])

def StarkZ_173_aBJ(L0,Sigma0,Omega0,J0,F10,F0,M0,L1,Sigma1,Omega1,J1,F11,F1,M1,S=1/2,I=5/2,iH=1/2):
    if not kronecker(L0,L1)*kronecker(Sigma0,Sigma1):
        return 0
    else:
        return (-1)**(F0-M0)*wigner_3j(F0,1,F1,-M0,0,M1)*\
            (-1)**(F1+F10+iH+1)*np.sqrt((2*F0+1)*(2*F1+1))*wigner_6j(F11,F1,iH,F0,F10,1)*\
            (-1)**(F11+J0+I+1)*np.sqrt((2*F10+1)*(2*F11+1))*wigner_6j(J1,F11,I,F10,J0,1)*\
            (-1)**(J0-Omega0)*wigner_3j(J0,1,J1,-Omega0,0,Omega1)*np.sqrt((2*J0+1)*(2*J1+1))


def TransitionDipole_173_aBJ(L0,Sigma0,Omega0,J0,F10,F0,M0,L1,Sigma1,Omega1,J1,F11,F1,M1,S=1/2,I=5/2,iH=1/2):
    if not kronecker(Sigma0,Sigma1):
        return 0
    else:
        TDM_total = sum([(-1)**(p)*\
            sum([
                (-1)**(F0-M0)*wigner_3j(F0,1,F1,-M0,p,M1)*(-1)**(F1+F10+iH+1)*np.sqrt((2*F0+1)*(2*F1+1))*\
                    wigner_6j(F11,F1,iH,F0,F10,1)*(-1)**(F11+J0+I+1)*np.sqrt((2*F10+1)*(2*F11+1))*\
                    wigner_6j(J1,F11,I,F10,J0,1)*(-1)**(J0-Omega0)*np.sqrt((2*J0+1)*(2*J1+1))*wigner_3j(J0,1,J1,-Omega0,q,Omega1)
                    for q in [-1,1]]) # Since L is changing, can only get q=+-1 transitions
            for p in range(-1,2)])
        # TDM_plus = sum([(-1)**(F0-M0)*wigner_3j(F0,1,F1,-M0,1,M1)*(-1)**(F1+J+I+1)*np.sqrt((2*F0+1)*(2*F1+1))*\
        #     wigner_6j(J1,F1,I,F0,J0,1)*(-1)**(J0-Omega0)*np.sqrt((2*J0+1)*(2*J1+1))*wigner_3j(J0,1,J1,-Omega0,q,Omega1) for q in range(-1,2)])
        # TDM_zero = sum([(-1)**(F0-M0)*wigner_3j(F0,1,F1,-M0,0,M1)*(-1)**(F1+J+I+1)*np.sqrt((2*F0+1)*(2*F1+1))*\
        #     wigner_6j(J1,F1,I,F0,J0,1)*(-1)**(J0-Omega0)*np.sqrt((2*J0+1)*(2*J1+1))*wigner_3j(J0,1,J1,-Omega0,q,Omega1) for q in range(-1,2)])
        # TDM_minus = sum([(-1)**(F0-M0)*wigner_3j(F0,1,F1,-M0,-1,M1)*(-1)**(F1+J+I+1)*np.sqrt((2*F0+1)*(2*F1+1))*\
        #     wigner_6j(J1,F1,I,F0,J0,1)*(-1)**(J0-Omega0)*np.sqrt((2*J0+1)*(2*J1+1))*wigner_3j(J0,1,J1,-Omega0,q,Omega1) for q in range(-1,2)])
        # TDM_total = TDM_plus + TDM_zero + TDM_minus
        return TDM_total

def TransitionDipole_173_aBJ_noM(L0,Sigma0,Omega0,J0,F10,F0,M0,L1,Sigma1,Omega1,J1,F11,F1,M1,S=1/2,I=5/2,iH=1/2):
    if not kronecker(Sigma0,Sigma1)*(not kronecker(L0,L1)):
        return 0
    elif not (kronecker(F0,F1) or kronecker(F0+1,F1) or kronecker(F0-1,F1)):
        return 0
    else:
        TDM_total = 1/np.sqrt(2*F1+1)*sum([(-1)**(F1+F10+iH+1)*np.sqrt((2*F0+1)*(2*F1+1))*\
                        wigner_6j(F11,F1,iH,F0,F10,1)*(-1)**(F11+J0+I+1)*np.sqrt((2*F10+1)*(2*F11+1))*\
                        wigner_6j(J1,F11,I,F10,J0,1)*(-1)**(J0-Omega0)*np.sqrt((2*J0+1)*(2*J1+1))*wigner_3j(J0,1,J1,-Omega0,q,Omega1)\
                        for q in [-1,1]]) # Since L is changing, can only get q=+-1 transitions
        # TDM_plus = sum([(-1)**(F0-M0)*wigner_3j(F0,1,F1,-M0,1,M1)*(-1)**(F1+J+I+1)*np.sqrt((2*F0+1)*(2*F1+1))*\
        #     wigner_6j(J1,F1,I,F0,J0,1)*(-1)**(J0-Omega0)*np.sqrt((2*J0+1)*(2*J1+1))*wigner_3j(J0,1,J1,-Omega0,q,Omega1) for q in range(-1,2)])
        # TDM_zero = sum([(-1)**(F0-M0)*wigner_3j(F0,1,F1,-M0,0,M1)*(-1)**(F1+J+I+1)*np.sqrt((2*F0+1)*(2*F1+1))*\
        #     wigner_6j(J1,F1,I,F0,J0,1)*(-1)**(J0-Omega0)*np.sqrt((2*J0+1)*(2*J1+1))*wigner_3j(J0,1,J1,-Omega0,q,Omega1) for q in range(-1,2)])
        # TDM_minus = sum([(-1)**(F0-M0)*wigner_3j(F0,1,F1,-M0,-1,M1)*(-1)**(F1+J+I+1)*np.sqrt((2*F0+1)*(2*F1+1))*\
        #     wigner_6j(J1,F1,I,F0,J0,1)*(-1)**(J0-Omega0)*np.sqrt((2*J0+1)*(2*J1+1))*wigner_3j(J0,1,J1,-Omega0,q,Omega1) for q in range(-1,2)])
        # TDM_total = TDM_plus + TDM_zero + TDM_minus
        return TDM_total
