import numpy as np
import sympy as sy
from functools import partial
from YbOH_parameters import params_general
from YbOH_matrix_elements import MQM_bBS,EDM_bBS,Sz_bBJ,T2QYb_bBS,b2a_matrix_174,decouple_b_174

def H_174X(q_numbers,params,matrix_elements,symbolic=True,E=0,B=0,M_values='all',precision=5):
    q_str = list(q_numbers)     # Get keys for quantum number dict
    if symbolic:
        Ez,Bz = sy.symbols('Ez Bz')
        size = len(q_numbers[q_str[0]])
        # Need to construct empty matrices to fill with matrix elements
        # Sympy does not like numpy arrays, so convert to list
        H = np.zeros((size,size)).tolist()
        #Iz = np.zeros((size,size)).tolist()
        #Sz = np.zeros((size,size)).tolist()
        for i in range(size):
            for j in range(size):
                # State out is LHS of bra ket, state in is RHS
                state_out = {q+'0':q_numbers[q][i] for q in q_str}
                state_in = {q+'1':q_numbers[q][j] for q in q_str}
                q_args = {**state_out,**state_in}
                elements = {term: sy.N(element(**q_args)) for term, element in matrix_elements.items()}
                # The Hamiltonian
                H[i][j] = params['Be']*elements['N^2'] + params['Gamma_SR']*elements['N.S'] + \
                    params['bF']*elements['I.S'] + params['c']/3*np.sqrt(6)*elements['T2_0(I,S)']
                if M_values!='none':
                    H[i][j]+=params['g_S']*params['mu_B']*elements['ZeemanZ']*Bz - params['muE_X']*elements['StarkZ']*Ez
                if params.get('q_lD') is not None:
                    H[i][j] += params['q_lD']/2*elements['l-doubling']
                # H[i][j] = round(H[i][j],precision)
                #Iz[i][j] = params['c']*elements['Iz']
                #Sz[i][j] = elements['Sz']
        # Need to construct IzSz term and add to Hamiltonian
        #H=matadd(H,matmult(Iz,Sz))
        # Create symbolic object
        H_symbolic = sy.Matrix(H)
        # Use symbolic object to create function that given E and B values, returns a numpy array
        H_func = sy.lambdify((Ez,Bz), H_symbolic, modules='numpy')
        return H_func,H_symbolic
    # Same as above, but fully numeric
    # else:
    #     Ez,Bz = [E,B]
    #     size = len(q_numbers[q_str[0]])
    #     H = np.zeros((size,size))
    #     Iz = np.zeros((size,size))
    #     Sz = np.zeros((size,size))
    #     for i in range(size):
    #         for j in range(size):
    #             state_out = {q+'0':q_numbers[q][i] for q in q_str}
    #             state_in = {q+'1':q_numbers[q][j] for q in q_str}
    #             q_args = {**state_out,**state_in}
    #             elements = {term: element(**q_args) for term, element in matrix_elements.items()}
    #             H[i,j] = params['Be']*elements['N^2'] + params['Gamma_SR']*elements['N.S'] + \
    #                 params['b']*elements['I.S'] + \
    #                 bend*params['q_lD']/2*lD_bBJ(*state_out,*state_in)+\
    #                 params['g_S']*params['mu_B']*Bz*elements['ZeemanZ'] - params['muE_X']*Ez*elements['StarkZ']
    #             if params.get('q_lD') is not None:
    #                 H[i,j] += params['q_lD']/2*elements['l-doubling']
    #             Iz[i,j] = elements['Iz']
    #             Sz[i,j] = elements['Sz']
    #     H = H + params['c']*(Iz@Sz)
        return H

# See documentation for H_174X
def H_174A(q_numbers,params,matrix_elements,symbolic=True,E=0,B=0,SO=0,M_values='all',precision=5):
    q_str = list(q_numbers)
    if symbolic:
        Ez,Bz = sy.symbols('Ez Bz')
        size = len(q_numbers[q_str[0]])
        H = np.zeros((size,size)).tolist()
        for i in range(size):
            for j in range(size):
                state_out = {q+'0':q_numbers[q][i] for q in q_str}
                state_in = {q+'1':q_numbers[q][j] for q in q_str}
                q_args = {**state_out,**state_in}
                elements = {term: element(**q_args) for term, element in matrix_elements.items()}
                H[i][j] = params['Be']*elements['N^2'] + SO*params['ASO']*elements['SO'] + \
                    params['bF']*elements['I.S'] + params['c']*np.sqrt(6)/3*elements['T2_0(IS)']+\
                    params['p+2q']*elements['Lambda-Doubling']
                if M_values!='none':
                    H[i][j]+=params['g_L']*params['mu_B']*Bz*elements['ZeemanLZ']+params['g_S']*params['mu_B']*Bz*elements['ZeemanSZ'] +\
                    Bz*params['g_lp']*params['mu_B']*elements['ZeemanParityZ'] - params['muE_A']*Ez*elements['StarkZ']
                # H[i][j] = round(H[i][j],precision)
        H_symbolic = sy.Matrix(H)
        H_func = sy.lambdify((Ez,Bz), H_symbolic, modules='numpy')
        return H_func,H_symbolic
    # else:
    #     Ez,Bz = [E,B]
    #     size = len(q_numbers[q_str[0]])
    #     H = np.zeros((size,size))
    #     Iz = np.zeros((size,size))
    #     Sz = np.zeros((size,size))
    #     for i in range(size):
    #         for j in range(size):
    #             state_out = {q+'0':q_numbers[q][i] for q in q_str}
    #             state_in = {q+'1':q_numbers[q][j] for q in q_str}
    #             q_args = {**state_out,**state_in}
    #             elements = {term: element(**q_args) for term, element in matrix_elements.items()}
    #             H[i,j] = params['Be']*elements['N^2'] + SO*params['ASO']*elements['SO'] + \
    #                 (params['bF']-params['c']/3)*elements['I.S'] + params['c']*elements['IzSz']+\
    #                 params['g_L']*params['mu_B']*Bz*elements['ZeemanLZ']+params['g_S']*params['mu_B']*Bz*elements['ZeemanSZ'] +\
    #                 Bz*params['g_lp']*params['mu_B']*elements['ZeemanParityZ'] - params['muE_A']*Ez*elements['StarkZ']+\
    #                 params['p+2q']*elements['Lambda-Doubling']
        return H

# See documentation for H_174X
def H_173X(q_numbers,params,matrix_elements,symbolic=True,E=0,B=0,M_values='all',precision=5):
    q_str = list(q_numbers)
    if symbolic:
        Ez,Bz = sy.symbols('Ez Bz')
        size = len(q_numbers[q_str[0]])
        H = np.zeros((size,size)).tolist()
        for i in range(size):
            for j in range(size):
                state_out = {q+'0':q_numbers[q][i] for q in q_str}
                state_in = {q+'1':q_numbers[q][j] for q in q_str}
                q_args = {**state_out,**state_in}
                elements = {term: element(**q_args) for term, element in matrix_elements.items()}
                H[i][j] = params['Be']*elements['N^2'] + params['Gamma_SR']*elements['N.S'] + \
                    params['bFYb']*elements['IYb.S'] + np.sqrt(6)*params['cYb']*elements['T2_0(IYb,S)'] +\
                    np.sqrt(6)/(4*5/2*(2*5/2-1))*params['e2Qq0']*elements['T2_0(IYb^2)'] +\
                    params['bFH']*elements['IH.S'] + (-np.sqrt(10))*params['cH']/3*elements['T2_0(IH,S)']
                if M_values != 'none':
                    H[i][j]+=params['g_S']*params['mu_B']*Bz*elements['ZeemanZ'] - params['muE_X']*Ez*elements['StarkZ']
                if params.get('q_lD') is not None:
                    H[i][j] += params['q_lD']/2*elements['l-doubling']
                # H[i][j] = round(H[i][j],precision)
        H_symbolic = sy.Matrix(H)
        H_func = sy.lambdify((Ez,Bz), H_symbolic, modules='numpy')
        return H_func,H_symbolic
    # else:
    #     Ez,Bz = [E,B]
    #     size = len(q_numbers[q_str[0]])
    #     H = np.zeros((size,size))
    #     for i in range(size):
    #         for j in range(size):
    #             state_out = {q+'0':q_numbers[q][i] for q in q_str}
    #             state_in = {q+'1':q_numbers[q][j] for q in q_str}
    #             q_args = {**state_out,**state_in}
    #             elements = {term: element(**q_args) for term, element in matrix_elements.items()}
    #             H[i,j] = params['Be']*elements['N^2'] + params['Gamma_SR']*elements['N.S'] + \
    #                 params['bFYb']*elements['IYb.S'] + np.sqrt(6)*params['cYb']*elements['T2_0(IYb,S)'] +\
    #                 np.sqrt(6)/(4*5/2*(2*5/2-1))*params['e2Qq0']*elements['T2_0(IYb^2)'] +\
    #                 params['bFH']*elements['IH.S'] + (-np.sqrt(10))*params['cH']/3*elements['T2_0(IH,S)']+\
    #                 params['g_S']*params['mu_B']*Bz*elements['ZeemanZ'] - params['muE_X']*Ez*elements['StarkZ']
    #             if params.get('q_lD') is not None:
    #                 H[i,j] += params['q_lD']/2*elements['l-doubling']
        return H

def H_173A(q_numbers,params,matrix_elements,symbolic=True,E=0,B=0,M_values='all',precision=5):
    q_str = list(q_numbers)
    if symbolic:
        Ez,Bz = sy.symbols('Ez Bz')
        size = len(q_numbers[q_str[0]])
        H = np.zeros((size,size)).tolist()
        for i in range(size):
            for j in range(size):
                state_out = {q+'0':q_numbers[q][i] for q in q_str}
                state_in = {q+'1':q_numbers[q][j] for q in q_str}
                q_args = {**state_out,**state_in}
                elements = {term: element(**q_args) for term, element in matrix_elements.items()}
                H[i][j] = params['Be']*elements['N^2'] + SO*params['ASO']*elements['SO']+\
                    params['h1/2Yb']*elements['IzLz_Yb'] + params['dYb']*elements['T2_2(IS)_Yb']+params['e2Qq0']*elements['T2_0(II)_Yb']+\
                    params['p+2q']*elements['Lambda-Doubling']
                if M_values!='none':
                    H[i][j]+=params['g_L']*params['mu_B']*Bz*elements['ZeemanLZ']+params['g_S']*params['mu_B']*Bz*elements['ZeemanSZ'] +\
                    Bz*params['g_lp']*params['mu_B']*elements['ZeemanParityZ'] - params['muE_A']*Ez*elements['StarkZ']
                # H[i][j] = round(H[i][j],precision)


                    # params['bFH']*elements['I.S'] + params['cH']*np.sqrt(6)/3*elements['T2_0(IS)']
        H_symbolic = sy.Matrix(H)
        H_func = sy.lambdify((Ez,Bz), H_symbolic, modules='numpy')
        return H_func,H_symbolic
    else:
        Ez,Bz = [E,B]
        size = len(q_numbers[q_str[0]])
        H = np.zeros((size,size))
        Iz = np.zeros((size,size))
        Sz = np.zeros((size,size))
        for i in range(size):
            for j in range(size):
                state_out = {q+'0':q_numbers[q][i] for q in q_str}
                state_in = {q+'1':q_numbers[q][j] for q in q_str}
                q_args = {**state_out,**state_in}
                elements = {term: element(**q_args) for term, element in matrix_elements.items()}
                H[i,j] = params['Be']*elements['N^2'] + SO*params['ASO']*elements['SO'] #+ \
                    # (params['bF']-params['c']/3)*elements['I.S'] + params['c']*elements['IzSz']+\
                    # params['g_L']*params['mu_B']*Bz*elements['ZeemanLZ']+params['g_S']*params['mu_B']*Bz*elements['ZeemanSZ'] +\
                    # Bz*params['g_lp']*params['mu_B']*elements['ZeemanParityZ'] - params['muE_A']*Ez*elements['StarkZ']+\
                    # params['p+2q']*elements['Lambda-Doubling']
        return H

def build_PTV_bBS(q_numbers,EDM_or_MQM):
    q_str = list(q_numbers)
    size = len(q_numbers[q_str[0]])
    H_PTV = np.zeros((size,size))
    EDM = {False: 0, True: 1}[EDM_or_MQM=='EDM']
    Mzz = {False: 0, True: 1}[EDM_or_MQM=='MQM']
    for i in range(size):
        for j in range(size):
            state_out = {q+'0':q_numbers[q][i] for q in q_str}
            state_in = {q+'1':q_numbers[q][j] for q in q_str}
            q_args = {**state_out,**state_in}
            H_PTV[i,j] = -Mzz/(2*5/2*(2*5/2-1))*(-2/3*np.sqrt(15))*MQM_bBS(**q_args)-EDM*EDM_bBS(**q_args)
    #         elif H_MQM==2:
    #             H1[i,j] = -Mzz/(2*5/2*(2*5/2-1))*2/3*np.sqrt(6)*T2QYb_bBS(**q_args)
    #             H2[i,j] = EDM_bBS(**q_args)
    #             H_PTV[i,j] = -EDM*EDM_bBS(**q_args)
    # if H_MQM==2:
    #     H_PTV = H_PTV + H1@H2
    return H_PTV

def build_PTV_bBJ(q_numbers):
    q_str= list(q_numbers)
    size = len(q_numbers[q_str[0]])
    H_PTV = np.zeros((size,size))
    EDM = 1
    for i in range(size):
        for j in range(size):
            state_out = {q+'0':q_numbers[q][i] for q in q_str}
            state_in = {q+'1':q_numbers[q][j] for q in q_str}
            q_args = {**state_out,**state_in}
            H_PTV[i,j] = -EDM*Sz_bBJ(**q_args)
    return H_PTV

def build_TDM_aBJ(q_in,q_out,TDM_matrix_element):
    q_str_in = list(q_in)
    q_str_out = list(q_out)
    size_in = len(q_in[q_str_in[0]])
    size_out = len(q_out[q_str_out[0]])
    H_TDM = np.zeros((size_out,size_in))
    for i in range(size_out):
        for j in range(size_in):
            state_out = {q+'0':q_out[q][i] for q in q_out}
            state_in = {q+'1':q_in[q][j] for q in q_in}
            q_args = {**state_out,**state_in}
            H_TDM[i,j] = TDM_matrix_element(**q_args)
    return H_TDM

def matmult(a,b):
    zip_b = list(zip(*b))
    return [[sum(ele_a*ele_b for ele_a, ele_b in zip(row_a, col_b))
             for col_b in zip_b] for row_a in a]

def matadd(a,b):
    return [[ele_a+ele_b for ele_a,ele_b in zip(row_a,row_b)] for row_a,row_b in zip(a,b)]

def convert_ab(input_qnumbers,output_qnumbers,S=1/2):
    input_keys = list(input_qnumbers)
    output_keys = list(output_qnumbers)
    input_size = len(input_qnumbers[input_keys[0]])
    output_size = len(output_qnumbers[output_keys[0]])
    basis_matrix = np.zeros((output_size,input_size))
    for i in range(output_size):
        for j in range(input_size):
            if 'N' in input_keys: #Convert case (b) to (a)
                a_qnumbers = {q:output_qnumbers[q][i] for q in output_keys}
                b_qnumbers = {q:input_qnumbers[q][j] for q in input_keys}
            else:
                b_qnumbers = {q:output_qnumbers[q][i] for q in output_keys}
                a_qnumbers = {q:input_qnumbers[q][j] for q in input_keys}
            basis_matrix[i,j] = b2a_matrix_174(a_qnumbers,b_qnumbers,S=S)
    return basis_matrix

def decouple_b(input_qnumbers,output_qnumbers,S=1/2):
    input_keys = list(input_qnumbers)
    output_keys = list(output_qnumbers)
    input_size = len(input_qnumbers[input_keys[0]])
    output_size = len(output_qnumbers[output_keys[0]])
    basis_matrix = np.zeros((output_size,input_size))
    for i in range(output_size):
        for j in range(input_size):
            decoupled_qnumbers = {q:output_qnumbers[q][i] for q in output_keys}
            b_qnumbers = {q:input_qnumbers[q][j] for q in input_keys}
            basis_matrix[i,j] = decouple_b_174(decoupled_qnumbers,b_qnumbers,S=S)
    return basis_matrix
