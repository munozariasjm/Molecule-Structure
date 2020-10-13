from YbOH_Library_class import YbOH_Library
import numpy as np
import sympy as sy
from numpy import linalg as npLA
from scipy import linalg as sciLA
import matplotlib.pyplot as plt
from copy import deepcopy

# Initialize a library with relevant functions and parameters for all states
library = YbOH_Library()

class YbOHLevels(object):

    '''This class is used to determine energy levels for a given vibronic state
    in YbOH. As an input, the user must specify the isotope (string), the state
    (string), and the range of N values (tuple: (Nmin, Nmax)).

    Example isotope formatting: '174'
    Example state formatting: 'X000'

    All calculations are done in MHz, with E in V/cm, and B in Gauss
    '''

    @classmethod
    def initialize_state(cls,isotope,state,N_range,all_M=True,round=5):
        if isotope not in ['170','171','172','173','174','176']:
            print(isotope, ' is not a valid isotope of Yb')
            return None
        if isotope not in ['173','174']:
            print(isotope, ' is an isotope not yet supported by this code')
            return None
        if state not in ['X000','X010','A000']:
            print('Your input state, ', state, ' is not currently supported by this code. \nAn example state string: X000')
            return None
        iso_state = isotope + state

        # Properties contain information relevant to the isotope and state of interest
        properties = {
            'iso_state': iso_state,
            'isotope': isotope,
            'state': state,
            'parameters': library.parameters[iso_state], # Hamiltonian parameters relevant to state and isotope
            'matrix_elements': library.matrix_elements[iso_state], #
            'hunds_case': library.cases[iso_state],
            'Lambda': library.Lambda[iso_state],
            'N_range': N_range,
            'all_M': all_M,
            'round': round #how much to round eigenvalues and eigenvectors

            # For later implementation:
            #'e_spin': electron_spin[iso_state],
            #'IYb_spin': nuclear_spin[iso_state]['Yb'],
            #'iH_spin': nuclear_spin[iso_state]['H'],
        }
        return cls(**properties)


    def __init__(self, **properties):
        # Create attributes using properties dict
        self.__dict__.update(properties)

        # Create quantum number dictionary, contains arrays of angular momenta eigenvalues indexed by basis vector
        # Example: {'J':[0,1,1,1], 'M':[0,-1,0,1]}
        self.q_numbers = library.q_number_builders[self.iso_state](self.N_range, all_M = self.all_M)
        self.q_str = list(self.q_numbers)

        # Create quantum numbers for alternate bases, like decoupled basis
        self.alt_q_numbers = {basis: q_builder(self.N_range,all_M=self.all_M)
            for basis,q_builder in library.alt_q_number_builders[self.iso_state].items()}

        # Create Hamiltonian.
        # H_function is a function that takes (E,B) values as an argument, and returns a numpy matrix
        # H_symbolic is a symbolic sympy matrix of the Hamiltonian
        self.H_function, self.H_symbolic = library.H_builders[self.iso_state](self.q_numbers)

        # Find free field eigenvalues and eigenvectors
        self.eigensystem(0,1e-4)

        # These attrbiutes will be used to store Zeeman and Stark information
        # Each index corresponds to a value of E or B.
        # Evals contains an array of eigenvalues for each field value
        # Evecs contains an array of eigenvectors for each field value
        self.Ez = None
        self._Bz = None
        self.evals_E = None
        self.evecs_E = None

        self.Bz = None
        self._Ez = None
        self.evals_B = None
        self.evecs_B = None

        # Attributes used for evaluating PT violating shifts
        self.H_PTV = None
        self.PTV_E = None
        self.PTV0 = None
        self.PTV_type = None


    def eigensystem(self,Ez_val,Bz_val,method='scipy_h',order=True, set_attr=True):
        if method == 'numpy':
            w,v = npLA.eig(0.5*(self.H_function(Ez_val, Bz_val)+np.transpose(self.H_function(Ez_val, Bz_val))))
        elif method == 'numpy_h':
            w,v = npLA.eigh(self.H_function(Ez_val, Bz_val))
        elif method == 'scipy':
            w,v = sciLA.eig(0.5*(self.H_function(Ez_val, Bz_val)+np.transpose(self.H_function(Ez_val, Bz_val))))
        elif method == 'scipy_h':
            w,v, = sciLA.eigh(self.H_function(Ez_val, Bz_val))
        evals = np.real(w)
        evecs = np.array([np.round(np.real(v[:,i]),self.round) for i in range(len(w))])
        for i,evec in enumerate(evecs):
            evecs[i]/=evec@evec
        if order:
            evals,evecs = order_eig(evals,evecs)
        if set_attr:
            self.evals0,self.evecs0 = [evals,evecs]
            self.E0,self.B0 = [Ez_val, Bz_val]
        return evals,evecs

    # Bz must be in Gauss
    def ZeemanMap(self,Bz_array,Ez_val=0,plot=False,**kwargs):
        self.Bz = Bz_array
        self._Ez = Ez_val
        evals0,evecs0 =self.eigensystem(Ez_val,Bz_array[0],set_attr=False)
        evals_B = [evals0]
        evecs_B = [evecs0]
        evecs_old = evecs0
        for B in Bz_array[1:]:
            evals_new,evecs_new = self.eigensystem(Ez_val,B)
            order = state_ordering(evecs_old,evecs_new,round=self.round-1)
            evecs_ordered = evecs_new[order,:]
            evals_ordered = evals_new[order]
            evecs_B.append(evecs_ordered)
            evals_B.append(evals_ordered)
            evecs_old = evecs_B[-1]
        evals_B,evecs_B = [np.array(evals_B),np.array(evecs_B)]
        self.evals_B = evals_B
        self.evecs_B = evecs_B
        if plot:
            self.plot_evals_EB('B',**kwargs)
        return

    # Ez must be in V/cm
    def StarkMap(self,Ez_array,Bz_val=1e-9,plot=False,**kwargs):
        self.Ez = Ez_array
        self._Bz = Bz_val
        evals0,evecs0 = self.eigensystem(Ez_array[0],Bz_val,set_attr=False)
        evals_E = [evals0]
        evecs_E = [evecs0]
        evecs_old = evecs0
        for E in Ez_array[1:]:
            evals_new,evecs_new = self.eigensystem(E,Bz_val)
            order = state_ordering(evecs_old,evecs_new,round=self.round-1)
            evecs_ordered = evecs_new[order,:]
            evals_ordered = evals_new[order]
            evecs_E.append(evecs_ordered)
            evals_E.append(evals_ordered)
            evecs_old = evecs_E[-1]
        evals_E,evecs_E = [np.array(evals_E),np.array(evecs_E)]
        self.evals_E = evals_E
        self.evecs_E = evecs_E
        if plot:
            self.plot_evals_EB('E',**kwargs)
        return

    def g_eff(self,step=1e-7):
        return g_eff_evecs(self.evals,self.evecs,self.E0,self.B0,step=step)

    def g_eff_evecs(self,evals,evecs,Ez,Bz,step=1e-7):
        evals0,evecs0 = evals,evecs
        evals1,evecs1 = self.eigensystem(Ez,Bz+step,set_attr=False)
        order = state_ordering(evecs0,evecs1,round=self.round-1)
        evecs1_ordered = evecs1[order,:]
        evals1_ordered = evals1[order]
        g_eff = []
        for E0,E1 in zip(evals0,evals1_ordered):
            g_eff.append((E1-E0)/(step*self.parameters['mu_B']))
        g_eff = np.array(g_eff)
        return g_eff


    def g_eff_EB(self,Ez,Bz,step=1e-7):
        self.eigensystem(Ez,Bz)
        g_eff = self.g_eff(step=step)
        return g_eff

    def PTV_shift(self,EDM_or_MQM):
        if '174' in self.iso_state:
            self.PTV_type = 'EDM'
            H_PTV = library.PTV_builders[self.iso_state](self.q_numbers)
        elif '173' in self.iso_state:
            self.PTV_type = EDM_or_MQM
            H_PTV = library.PTV_builders[self.iso_state](self.q_numbers, EDM_or_MQM)
        self.H_PTV = H_PTV
        evals,evecs = self.evals0,self.evecs0
        PTV_shift = []
        for evec in evecs:
            E_PTV = evec@H_PTV@evec
            PTV_shift.append(E_PTV)
        PTV_shift = np.array(PTV_shift)
        return PTV_shift

    def PTV_shift_EB(self,Ez,Bz,EDM_or_MQM):
        self.eigensystem(Ez,Bz)
        PTV_shift = self.PTV_shift(EDM_or_MQM)
        return PTV_shift

    def plot_evals_EB(self,E_or_B,kV_kG=False, GHz=False):
        x_scale = {False: 1, True: 10**-3}[kV_kG]
        y_scale = {False: 1, True: 10**-3}[GHz]

        field,evals = {
            'E': (self.Ez,self.evals_E),
            'B': (self.Bz,self.evals_B)
        }[E_or_B]

        evals = evals.T

        x_label = {
            True: {'E': 'E (kV/cm)', 'B': 'B (kGauss)'}[E_or_B],
            False: {'E': 'E (V/cm)', 'B': 'B (Gauss)'}[E_or_B]
        }[kV_kG]

        y_label = {
            True: 'Energy (GHz)',
            False: 'Energy (MHz)'
        }[GHz]

        state_str = {
            '174X000': r'$^{174}$YbOH $\tilde{X}(000)$',
            '174X010': r'$^{174}$YbOH $\tilde{X}(010)$',
            '173X000': r'$^{173}$YbOH $\tilde{X}(000)$',
            '173X010': r'$^{173}$YbOH $\tilde{X}(010)$',
            '174A000': r'$^{174}$YbOH $\tilde{A}(000)$',
            '173A000': r'$^{173}$YbOH $\tilde{A}(000)$'
        }[self.iso_state]

        EB_str = {
            'E': 'Stark Shifts',
            'B': 'Zeeman Shifts'
        }[E_or_B]

        title = state_str + ' ' + EB_str + r', $N={}$'.format(str(self.N_range)[1:-1])

        plt.figure(figsize=(10,7))
        for trace in evals:
            plt.plot(x_scale*field,y_scale*trace)
        plt.xlabel(x_label,fontsize=14)
        plt.ylabel(y_label,fontsize=14)
        plt.title(title,fontsize=16)
        return


    def write_state(self,eval_i,Ez=None,Bz=None):
        if Ez is None and Bz is None:
            pass
        else:
            if Ez==self.E0 and Bz==self.B0:
                pass
            else:
                self.eigensystem(Ez,Bz,set_attr=True)
        i=eval_i
        if i<0:
            i = len(self.evals0)+i
        vector=self.evecs0[i]
        energy = self.evals0[i]
        print('E = {} MHz\n'.format(energy))
        if self.PTV0 is not None:
            print('{} Shift = {}\n'.format(self.PTV_type,self.PTV0[i]))
        #sum_sq = 0
        for index in np.nonzero(vector)[0]:
            v={q:self.q_numbers[q][index] for q in self.q_numbers}
            coeff = vector[index]
            if self.hunds_case == 'bBS':
                print(' {} |\u039B={},N={},G={},F1={},F={},M={}> \n'.format(coeff,v['L'],v['N'],v['G'],v['F1'],v['F'],v['M']))
            elif self.hunds_case == 'bBJ':
                print(' {} |\u039B={},N={},J={},F={},M={}> \n'.format(coeff,v['L'],v['N'],v['J'],v['F'],v['M']))
            elif self.hunds_case == 'aBJ':
                print(' {} |\u039B={},\u03A3={},\u03A9={},J={},F={},M={}> \n'.format(coeff,v['L'],v['Sigma'],v['Omega'],v['J'],v['F'],v['M']))

    def PTV_Map(self,EDM_or_MQM,E_or_B='E',plot=False):
        if '174' in self.iso_state:
            self.PTV_type = 'EDM'
            H_PTV = library.PTV_builders[self.iso_state](self.q_numbers)
        elif '173' in self.iso_state:
            self.PTV_type = EDM_or_MQM
            H_PTV = library.PTV_builders[self.iso_state](self.q_numbers, EDM_or_MQM)
        self.H_PTV = H_PTV
        if E_or_B=='E':
            PTV_vs_E = []
            if self.evecs_E is None:
                print('Run StarkMap first')
                return None
            for evecs in self.evecs_E:
                PTV_shift = []
                for evec0 in evecs:
                    E_PTV = evec0@H_PTV@evec0
                    PTV_shift.append(np.round(E_PTV,self.round))
                PTV_shift = np.array(PTV_shift)
                PTV_vs_E.append(PTV_shift)
            PTV_vs_E = np.array(PTV_vs_E)
            self.PTV_E = PTV_vs_E
        elif E_or_B=='B':
            PTV_vs_B = []
            if self.evecs_B is None:
                print('Run ZeemanMap first')
                return None
            for evecs in self.evecs_B:
                PTV_shift = []
                for evec0 in evecs:
                    E_PTV = evec0@H_PTV@evec0
                    PTV_shift.append(np.round(E_PTV,self.round))
                PTV_shift = np.array(PTV_shift)
                PTV_vs_B.append(PTV_shift)
            PTV_vs_B = np.array(PTV_vs_B)
            self.PTV_B = PTV_vs_B
        if plot:
            self.plot_PTV(E_or_B)
        return


    def g_eff_Map(self,E_or_B='E',step=1e-7):
        if E_or_B=='E':
            if self.evecs_E is None:
                print('Run StarkMap first')
                return None
            g_eff_vs_E = []
            for i,evecs in enumerate(self.evecs_E):
                g_eff = self.g_eff_evecs(self.evals_E[i],evecs,self.Ez[i],self._Bz,step=step)
                g_eff_vs_E.append(g_eff)
            g_eff_vs_E = np.array(g_eff_vs_E)
            self.g_eff_E = g_eff_vs_E
            return g_eff_vs_E
        else:
            if self.evecs_B is None:
                print('Run ZeemanMap first')
                return None
            g_eff_vs_B = []
            for i,evecs in enumerate(self.evecs_B):
                g_eff = self.g_eff_evecs(self.evals_B[i],evecs,self._Ez,self.Bz[i],step=step)
                g_eff_vs_B.append(g_eff)
            g_eff_vs_B = np.array(g_eff_vs_B)
            self.g_eff_B = g_eff_vs_B
            return g_eff_vs_B


    def plot_PTV(self,E_or_B='E',kV_kG=False):

        if self.PTV_E is None and E_or_B=='E':
            print('Need to run PTV_Map first')
            return
        if self.PTV_B is None and E_or_B == 'B':
            print('Need to run PTV_BMap first')
            return

        x_scale = {False: 1, True: 10**-3}[kV_kG]

        x_label = {
            True: {'E': 'E (kV/cm)', 'B': 'B (kGauss)'}[E_or_B],
            False: {'E': 'E (V/cm)', 'B': 'B (Gauss)'}[E_or_B]
        }[kV_kG]

        field,shifts = {'E': [self.Ez,self.PTV_E], 'B':[self.Bz,self.PTV_B]}[E_or_B]
        shifts = shifts.T # change primary index from E field to eigenvector

        y_label = 'Shift Energy (MHz)'

        state_str = {
            '174X000': r'$^{174}$YbOH $\tilde{X}(000)$',
            '174X010': r'$^{174}$YbOH $\tilde{X}(010)$',
            '173X000': r'$^{173}$YbOH $\tilde{X}(000)$',
            '173X010': r'$^{173}$YbOH $\tilde{X}(010)$',
            '174A000': r'$^{174}$YbOH $\tilde{A}(000)$',
            '173A000': r'$^{173}$YbOH $\tilde{A}(000)$'
        }[self.iso_state]

        PTV_str = {
            'EDM': 'EDM Shifts',
            'MQM': 'MQM Shifts'
        }[self.PTV_type]


        title = state_str + ' ' + PTV_str + r', $N={}$'.format(str(self.N_range)[1:-1])

        plt.figure(figsize=(10,7))
        for trace in shifts:
            plt.plot(x_scale*field,trace)
        plt.xlabel(x_label,fontsize=14)
        plt.ylabel(y_label,fontsize=14)
        plt.title(title,fontsize=16)
        return

    def filter_evecs(self,evecs,q_filter,filter_vals):
        idx_filtered = []
        for i,evec in enumerate(evecs):
            max_idx = np.argmax(evec**2)
            for val in filter_vals:
                if self.q_numbers[q_filter][max_idx]==val:
                    idx_filtered.append(i)
        idx_filtered=np.array(idx_filtered)
        return idx_filtered

    def display_levels(self,Ez,Bz,pattern_q,label=True,label_q =None,width=0.75,figsize=(10,10),ylim=None):
        if label_q is None:
            label_q = self.q_str
        if Ez==self.E0 and Bz==self.B0:
            evals,evecs = [self.evals0,self.evecs0]
        else:
            evals,evecs = self.eigensystem(Ez,Bz, set_attr=True)
        if ylim is None:
            scale = abs(evals[-1]-evals[0])
            ylim = (evals[0]-0.1*scale,evals[-1]+0.1*scale)
        else:
            scale = abs(ylim[1]-ylim[0])
        fig = plt.figure(figsize=figsize)
        primary_q = {q:[] for q in self.q_str}
        M_bounds = []
        for evec in evecs:
            for q in label_q:
                max_q_val = self.q_numbers[q][np.argmax(evec**2)]
                primary_q[q].append(max_q_val)
                if q=='M':
                    M_bounds.append([(max_q_val-width/2),(max_q_val+width/2)])
        M_bounds = np.array(M_bounds).T
        plt.hlines(evals,*M_bounds)
        plt.ylim(ylim)
    #     left,right = plt.xlim()
    #     plt.xlim(left-1,right+2)
        if label:
            off = 0.03
            sign = 1
            label_q_no_M = [q for q in label_q if q!='M']
            prev_pattern = None
            for i,energy in enumerate(evals):
                current_pattern = primary_q[pattern_q][i]
                if prev_pattern != current_pattern and (ylim[0]+scale*off < energy < ylim[1]-scale*off):
                    label_str = r'$|$'
                    for q in label_q_no_M:
                        label_str+=r'${}$={},'.format(q,primary_q[q][i])
                    label_str = label_str[:-1]+r'$\rangle$'
                    plt.annotate(label_str,(0,energy-sign*off*scale),ha='center',fontsize=10)
                    if Ez ==0: # alternating label location is only useful at zero field
                        sign*=-1
                prev_pattern = current_pattern
            bot,top = plt.ylim()
            plt.ylim(bot-0.03*scale,top+0.03*scale)
        plt.xlabel(r'$M_F$',fontsize=14)
        plt.ylabel('Energy (MHz)',fontsize=14)

        return

    def display_PTV(self,Ez,Bz,EDM_or_MQM,width=0.75,figsize=(9,9),ylim=None):
        if '174' in self.iso_state:
            self.PTV_type = 'EDM'
            H_PTV = library.PTV_builders[self.iso_state](self.q_numbers)
        elif '173' in self.iso_state:
            self.PTV_type = EDM_or_MQM
            H_PTV = library.PTV_builders[self.iso_state](self.q_numbers, EDM_or_MQM)
        if Ez==self.E0 and Bz==self.B0:
            evals,evecs = [self.evals0,self.evecs0]
        else:
            evals,evecs = self.eigensystem(Ez,Bz, set_attr=True)
        PTV_shifts = []
        for evec0 in evecs:
            E_PTV = evec0@H_PTV@evec0 #should I complex conjugate?
            PTV_shifts.append(np.round(E_PTV,self.round))
        PTV_shifts = np.array(PTV_shifts)
        self.PTV0 = PTV_shifts
        if ylim is None:
            scale = abs(evals[-1]-evals[0])
            ylim = (evals[0]-0.1*scale,evals[-1]+0.1*scale)
        else:
            scale = abs(ylim[1]-ylim[0])
        fig = plt.figure(figsize=figsize)
        plt.ylim(ylim)
        M_bounds = []
        M_vals = []
        for evec in evecs:
            M = self.q_numbers['M'][np.argmax(evec**2)]
            M_vals.append(M)
            M_bounds.append([(M-width/2),(M+width/2)])
        M_bounds = np.array(M_bounds).T
        plt.hlines(evals,*M_bounds)
        off = 0.03
        for i,shift in enumerate(PTV_shifts):
            if (ylim[0]+scale*off < evals[i] < ylim[1]-scale*off):
                plt.annotate(shift,(M_vals[i],evals[i]-off*scale),ha='center',fontsize=12)
        bot,top = plt.ylim()
        plt.ylim(bot-0.3,top+0.3)
        plt.xlabel(r'$M_F$',fontsize=14)
        plt.ylabel('Energy (MHz)',fontsize=14)
        return

    def convert_evecs(self,basis,evecs=None):
        if evecs is None:
            evecs = self.evecs0
        current_case = self.hunds_case
        new_case = basis
        if new_case == current_case:
            print('Eigenvectors already in {} basis'.format(current_case))
            return evecs
        inputt = self.q_numbers
        output = self.alt_q_numbers[new_case]
        if ('a' in new_case and 'b' in current_case) or ('b' in new_case and 'a' in current_case):
            basis_matrix = library.basis_changers['a_b'](inputt,output)
        elif ('decoupled' in new_case and 'b' in current_case):
            basis_matrix = library.basis_changers['b_decoupled'](inputt,output)
        elif ('decoupled' in new_case and 'a' in current_case):
            basis_matrix = library.basis_changers['a_decoupled'](inputt,output)
        converted_evecs = []
        for i in range(len(evecs)):
            converted_evecs.append(basis_matrix@evecs[i])
        converted_evecs = np.array(converted_evecs)
        print('Successfully converted eigenvectprs from {} to {}'.format(current_case,new_case))
        return converted_evecs

    def gen_state_str(self,vector,basis=None,label_q=None,thresh=0.01,show_coeff=True,new_line=False,round=None):
        q_numbers = self.q_numbers
        label_q = self.q_str
        if round is None:
            round=self.round
        if basis is not None:
            if 'decouple' in basis:
                q_numbers = self.alt_q_numbers['decoupled']
                label_q = list(self.alt_q_numbers['decoupled'])
            elif 'a' in basis:
                q_numbers = self.alt_q_numbers['aBJ']
                label_q = list(self.alt_q_numbers['aBJ'])
        full_label = r''
        nonzero_idx = np.nonzero(vector)[0]
        first = 0
        for i,index in enumerate(nonzero_idx):
            coeff = np.round(vector[index],round)
            if abs(coeff) < thresh:
                first+=1
                continue
            sign = {True: '+', False: '-'}[coeff > 0]
            sign0 = {True: ' ', False: sign}[sign=='+']
            if show_coeff:
                label_str = r'${}{}|'.format({True: sign0, False: sign}[i==first],abs(coeff))
            else:
                label_str = r'${}|'.format({True: sign0, False: sign}[i==first])
            if new_line:
                label_str = r'$'+label_str
            val = {q:q_numbers[q][index] for q in label_q}
            for q in label_q:
                _q = {True: '\u039B', False: q}[q=='L']
                _q = {True: '\u03A3', False: _q}[q=='Sigma']
                _q = {True: '\u03A9', False: _q}[q=='Omega']
                if (abs(val[q]) % 1) !=0:
                    label_str+=r'{}=\tfrac{{{}}}{{{}}},'.format(_q,*val[q].as_integer_ratio())
                else:
                    label_str+=r'{}={},'.format(_q,int(val[q]))
            label_str = label_str[:-1]+r'\rangle$'
            if new_line:
                label_str+=r'$'
            full_label+=label_str
        if i==first and show_coeff==False:
            if new_line:
                full_label = r'$$'+full_label[3:]
            else:
                full_label = r'$'+full_label[2:]
        return full_label

def XA_branching_ratios(X,A,Ez,Bz): # must be in case a
    A.eigensystem(Ez,Bz)
    X.eigensystem(Ez,Bz)
    X_evecs_a = X.convert_evecs('aBJ')
    TDM_matrix = library.TDM_builders[A.iso_state](A.q_numbers,X.alt_q_numbers['aBJ'])
    BR_matrix = (X_evecs_a@TDM_matrix@A.evecs0.T)**2
    return BR_matrix


def state_ordering(evecs_old,evecs_new,round=3):
    overlap = abs(np.round(evecs_old@evecs_new.T,round))     #Essentially a matrix of the fidelities: |<phi|psi>|
    #calculate trace distance
    for o in overlap:
        for _o in o:
            if (_o>1):
                print('OVERLAP BIGGER THAN 1', _o)
    trace_dist = np.sqrt(1-np.square(overlap))
    ordering = np.array([trace_dist[i,:].argmin() for i in range(len(evecs_old))])
    return ordering

def order_eig(evals,evecs):
    order = np.argsort(evals)
    evecs_ordered = np.array(evecs)[order,:]
    evals_ordered = np.array(evals)[order]
    return evals_ordered,evecs_ordered
