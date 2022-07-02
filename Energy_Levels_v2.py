from molecule_library_class import Molecule_Library
import numpy as np
import sympy as sy
from numpy import linalg as npLA
from scipy import linalg as sciLA
import matplotlib.pyplot as plt
from copy import deepcopy
# import multiprocessing as mp
from functools import partial
from time import perf_counter
import torch

class QState(objtect):

    '''This class is mean to represent a single "ket" or quantum state.
    Its attrbiutes are a basis and coefficients in that basis.
    Also energy.
    '''

    def __init__(self,basis,coeff_vector,q_numbers):
        self.basis = basis
        self.evec = coeff_vector
        self.q_numbers = q_numbers
        self.string = self.make_string()
        self.latex = self.make_latex()

        

class MoleculeLevels(object):

    '''This class is used to determine energy levels for a given vibronic state
    in YbOH. As an input, the user must specify the isotope (string), the state
    (string), and the range of N values (tuple: (Nmin, Nmax)).

    Example isotope formatting: '174'
    Example state formatting: 'X000'

    All calculations are done in MHz, with E in V/cm, and B in Gauss
    '''

    @classmethod
    def initialize_state(cls,molecule,isotope,state,N_range,M_values='all',I=[0,1/2],S=1/2,P_values=[],M_range=[],round=6):
        if molecule=='YbOH':
            if isotope not in ['170','171','172','173','174','176']:
                print(isotope, 'is not a valid isotope of Yb')
                return None
            if isotope not in ['173','174','171']:
                print(isotope, 'is an isotope not yet supported by this code')
                return None
            if state not in ['X000','X010','A000']:
                print('Your input state, ', state, ' is not currently supported by this code. \nAn example state string: X000')
                return None
        elif molecule=='CaOH':
            if isotope not in ['40','42','43','44','46']:
                print(isotope, 'is not a stable isotope of Ca')
                return None
            if isotope not in ['40']:
                print(isotope, 'is an isotope not yet supported by this code')
                return None
            if state not in ['X010','X000','A000','B000']:
                print('Your input state, ', state, ' is not currently supported by this code. \nAn example state string: X000')
                return None
        iso_state = isotope + state
        if P_values == []:
            print('No P values provided, using P=1/2 as default')
            P_values=[1/2]

            # Properties contain information relevant to the isotope and state of interest
        properties = {
            'molecule': molecule,
            'iso_state': iso_state,
            'isotope': isotope,
            'state': state,
            'N_range': N_range,
            'M_values': M_values,
            'round': round,     #how much to round eigenvalues and eigenvectors
            'e_spin': S,    #electronic spin number
            'I_spins': I,    #spin of nuclei, [I_Yb, I_H]. I=0 means ignore
            'M_range': M_range,
            'P_values': P_values
        }
        return cls(**properties)


    def __init__(self, **properties):
        # Create attributes using properties dict
        self.__dict__.update(properties)

        # Initialize a library with relevant functions and parameters for all states
        self.library = Molecule_Library(self.molecule,self.I_spins,self.M_values,self.P_values)
        self.parameters = self.library.parameters[self.iso_state] # Hamiltonian parameters relevant to state and isotope
        self.matrix_elements = self.library.matrix_elements[self.iso_state]
        self.hunds_case = self.library.cases[self.iso_state]
        self.K: self.library.K[self.iso_state]

        # Create quantum number dictionary, contains arrays of angular momenta eigenvalues indexed by basis vector
        # Example: {'J':[0,1,1,1], 'M':[0,-1,0,1]}
        self.q_numbers = self.library.q_number_builders[self.iso_state](self.N_range, I_list=self.I_spins,M_values = self.M_values,M_range=self.M_range)
        self.q_str = list(self.q_numbers)

        # Create quantum numbers for alternate bases, like decoupled basis
        self.alt_q_numbers = {basis: q_builder(self.N_range,M_values = self.M_values,M_range=self.M_range)
            for basis,q_builder in self.library.alt_q_number_builders[self.iso_state].items()}

        # Create Hamiltonian.
        # H_function is a function that takes (E,B) values as an argument, and returns a numpy matrix
        # H_symbolic is a symbolic sympy matrix of the Hamiltonian
        self.H_function, self.H_symbolic = self.library.H_builders[self.iso_state](self.q_numbers,M_values=self.M_values,precision=self.round)

        # Find free field eigenvalues and eigenvectors
        if self.M_values == 'all' or self.M_values == 'pos':
            self.eigensystem(0,1e-6)
        else:
            self.eigensystem(0,0)
        self.size = len(self.evecs0)
        self.Parity_mat = self.library.all_parity[self.iso_state](self.q_numbers,self.q_numbers)
        self.generate_parities(self.evecs0)

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
        self.PTV_B = None
        self.PTV0 = None
        self.PTV_type = None

        self.state_str =  r'$^{{{iso}}}${mol} $\tilde{{{state}}}({vib})$'.format(iso=self.isotope,mol=self.molecule,state = self.state[:1],vib=self.state[1:])


    def eigensystem(self,Ez_val,Bz_val,method='torch',order=True, set_attr=True, Normalize=False):
        evals,evecs = diagonalize(self.H_function(Ez_val,Bz_val),method=method,order=order, Normalize=Normalize,round=self.round)
        if set_attr:
            self.evals0,self.evecs0 = [evals,evecs]
            self.E0,self.B0 = [Ez_val, Bz_val]
        return evals,evecs


    # Bz must be in Gauss
    def ZeemanMap(self,Bz_array,Ez_val=0,plot=False,output=False,method='torch',**kwargs):
        self.Bz = Bz_array
        self._Ez = Ez_val
        # t0 = perf_counter()
        B_matrices = np.array([self.H_function(Ez_val,Bz_val) for Bz_val in Bz_array])
        #results = [self.eigensystem(Ez_val,Bz_val,set_attr=False) for Bz_val in Bz_array]
        # pool = mp.Pool(mp.cpu_count())
        # worker = partial(diagonalize,order=True, Normalize=False,round=self.round)
        # results = pool.starmap(worker,[(self.H_function(Ez_val,Bz_val),) for Bz_val in Bz_array.tolist()])
        # pool.close()
        # pool.join()
        # t1 = perf_counter()
        #evals_B, evecs_B = list(zip(*results))
        # evals_B = np.array(evals_B)
        # evecs_B = np.array(evecs_B)
        evals_B,evecs_B = diagonalize_batch(B_matrices,method=method)
        # t2 = perf_counter()
        evecs_old = evecs_B[0]
        step = abs(Bz_array[1] - Bz_array[0])/abs(Bz_array).max()
        for i in range(len(Bz_array)):
            if i == 0:
                continue
            evals_new,evecs_new = evals_B[i], evecs_B[i]
            order = state_ordering(evecs_old,evecs_new,round=self.round,step=step)
            evecs_ordered = evecs_new[order,:]
            evals_ordered = evals_new[order]
            evecs_B[i] = evecs_ordered
            evals_B[i] = evals_ordered
            evecs_old = evecs_B[i]
        # t3 = perf_counter()
        # print('1:',(t1-t0))
        # print('2:',(t2-t1))
        # print('3:',(t3-t2))
        # evals_B,evecs_B = [np.array(evals_B),np.array(evecs_B)]
        self.evals_B = evals_B
        self.evecs_B = evecs_B
        if plot:
            self.plot_evals_EB('B',**kwargs)
        if output:
            return evals_B,evecs_B
        else:
            return

    # Ez must be in V/cm
    def StarkMap(self,Ez_array,Bz_val=1e-9,plot=False,output=False,method='torch',**kwargs):
        self.Ez = Ez_array
        self._Bz = Bz_val
        E_matrices = np.array([self.H_function(Ez_val,Bz_val) for Ez_val in Ez_array])
        # results = [self.eigensystem(Ez_val,Bz_val,set_attr=False) for Ez_val in Ez_array.tolist()]
        # pool = mp.Pool(mp.cpu_count())
        # worker = partial(diagonalize,order=True, Normalize=False,round=self.round)
        # results = pool.starmap(worker,[(self.H_function(Ez_val,Bz_val),) for Ez_val in Ez_array.tolist()])
        # pool.close()
        # pool.join()
        # evals_E, evecs_E = list(zip(*results))
        # evals_E = np.array(evals_E)
        # evecs_E = np.array(evecs_E)
        evals_E,evecs_E = diagonalize_batch(E_matrices,method=method)
        evecs_old = evecs_E[0]
        step = abs(Ez_array[1] - Ez_array[0])/abs(Ez_array).max()
        for i in range(len(Ez_array)):
            if i == 0:
                continue
            evals_new,evecs_new = evals_E[i], evecs_E[i]
            order = state_ordering(evecs_old,evecs_new,round=self.round,step=step)
            evecs_ordered = evecs_new[order,:]
            evals_ordered = evals_new[order]
            evecs_E[i] = evecs_ordered
            evals_E[i] = evals_ordered
            evecs_old = evecs_E[i]
        evals_E,evecs_E = [np.array(evals_E),np.array(evecs_E)]
        self.evals_E = evals_E
        self.evecs_E = evecs_E
        if plot:
            self.plot_evals_EB('E',**kwargs)
        if output:
            return evals_E,evecs_E
        else:
            return

    def g_eff_EB(self,Ez=None,Bz=None,step=1e-7):
        recalc = False
        if Ez is None:
            Ez = self.E0
            recalc=True
        if Bz is None:
            Bz = self.B0
            recalc=True
        if recalc:
            self.eigensystem(Ez,Bz);
        return self.g_eff_evecs(self.evals0,self.evecs0,Ez,Bz,step=step)

    def g_eff_evecs(self,evals,evecs,Ez,Bz,step=1e-7):
        evals0,evecs0 = evals,evecs
        evals1,evecs1 = self.eigensystem(Ez,Bz+step,set_attr=False)
        order = state_ordering(evecs0,evecs1,round=self.round)
        evecs1_ordered = evecs1[order,:]
        evals1_ordered = evals1[order]
        g_eff = []
        for E0,E1 in zip(evals0,evals1_ordered):
            g_eff.append((E1-E0)/(step*self.parameters['mu_B']))
        g_eff = np.array(g_eff)
        return g_eff

#Is this redundant?
    # def g_eff_EB(self,Ez,Bz,step=1e-7):
    #     self.eigensystem(Ez,Bz)
    #     g_eff = self.g_eff(step=step)
    #     return g_eff

    def PTV_shift(self,EDM_or_MQM):
        if '174' in self.iso_state or '40' in self.iso_state:
            self.PTV_type = 'EDM'
            H_PTV = self.library.PTV_builders[self.iso_state](self.q_numbers)
        elif '173' in self.iso_state:
            self.PTV_type = EDM_or_MQM
            H_PTV = self.library.PTV_builders[self.iso_state](self.q_numbers, EDM_or_MQM)
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

    def plot_evals_EB(self,E_or_B,kV_kG=False, GHz=False,Freq=True):
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
        state_str = self.state_str

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


    def write_state(self,eval_i,Ez=None,Bz=None,show_PTV=False):
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
        if self.PTV0 is not None and show_PTV:
            print('{} Shift = {}\n'.format(self.PTV_type,self.PTV0[i]))
        #sum_sq = 0
        for index in np.nonzero(vector)[0]:
            v={q:self.q_numbers[q][index] for q in self.q_numbers}
            coeff = vector[index]
            if self.hunds_case == 'bBS':
                print(' {} |K={},N={},G={},F1={},F={},M={}> \n'.format(coeff,v['K'],v['N'],v['G'],v['F1'],v['F'],v['M']))
            elif self.hunds_case == 'bBJ':
                print(' {} |K={},N={},J={},F={},M={}> \n'.format(coeff,v['K'],v['N'],v['J'],v['F'],v['M']))
            elif self.hunds_case == 'aBJ':
                print(' {} |K={},\u03A3={},P={},J={},F={},M={}> \n'.format(coeff,v['K'],v['Sigma'],v['P'],v['J'],v['F'],v['M']))


    def PTV_Map(self,EDM_or_MQM,E_or_B='E',plot=False):
        if '174' in self.iso_state or '40' in self.iso_state:
            self.PTV_type = 'EDM'
            H_PTV = self.library.PTV_builders[self.iso_state](self.q_numbers)
        elif '173' in self.iso_state:
            self.PTV_type = EDM_or_MQM
            H_PTV = self.library.PTV_builders[self.iso_state](self.q_numbers, EDM_or_MQM)
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

        y_label = {
            'EDM': r'$\langle \Sigma \rangle$',
            'MQM': r'$\langle S \cdot T_{II} \cdot n \rangle'
        }[self.PTV_type]

        state_str =  self.state_str

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

    def display_levels(self,Ez,Bz,pattern_q,idx = None,label=True,parity=False,label_q =None,width=0.75,figsize=(10,10),ylim=None, deltaE_label = 3000,alt_label=False):
        if label_q is None:
            label_q = self.q_str
        if 'M' not in label_q:
            label_q.append('M')
        if pattern_q not in label_q:
            label_q.append(pattern_q)
        if Ez==self.E0 and Bz==self.B0:
            evals,evecs = [self.evals0,self.evecs0]
        else:
            evals,evecs = self.eigensystem(Ez,Bz, set_attr=False)
        if idx is not None:
            evals = evals[idx]
            evecs = evecs[idx]
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
            prev_energy = 0
            labeled_energies = []
            for i,energy in enumerate(evals):
                current_pattern = primary_q[pattern_q][i]
                if (prev_pattern != current_pattern or (prev_pattern == current_pattern and abs(energy-prev_energy)>deltaE_label)) and (ylim[0]+scale*off < energy < ylim[1]-scale*off):
                    if energy not in labeled_energies:
                        label_str = r'$|$'
                        if parity:
                            label_str+='${}$,'.format({1: '+', -1:'-'}[self.parities[i+idx[0]]])
                        for q in label_q_no_M:
                            label_str+=r'${}$={},'.format(q,primary_q[q][i])
                        label_str = label_str[:-1]+r'$\rangle$'
                        if self.M_values=='custom':
                            loc = primary_q['M'][i]
                        elif self.M_values=='positive':
                            loc = abs(primary_q['M'][i] % 1)
                        else:
                            loc = 0
                        plt.annotate(label_str,(loc,energy-sign*off*scale),ha='center',fontsize=10)
                        labeled_energies.append(energy)
                        if alt_label: # alternate labels
                            sign*=-1
                prev_energy = energy
                prev_pattern = current_pattern
            bot,top = plt.ylim()
            plt.ylim(bot-0.03*scale,top+0.03*scale)
        plt.xlabel(r'$M_F$',fontsize=14)
        plt.ylabel('Energy (MHz)',fontsize=14)

        return

    def display_PTV(self,Ez,Bz,EDM_or_MQM,width=0.75,figsize=(9,9),ylim=None,round=None):
        if '174' in self.iso_state or '40' in self.iso_state:
            self.PTV_type = 'EDM'
            H_PTV = self.library.PTV_builders[self.iso_state](self.q_numbers)
        elif '173' in self.iso_state:
            self.PTV_type = EDM_or_MQM
            H_PTV = self.library.PTV_builders[self.iso_state](self.q_numbers, EDM_or_MQM)
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
                if round is not None:
                    shift = np.round(shift,round)
                plt.annotate(shift,(M_vals[i],evals[i]-off*scale),ha='center',fontsize=12)
        bot,top = plt.ylim()
        plt.ylim(bot-0.3,top+0.3)
        plt.xlabel(r'$M_F$',fontsize=14)
        plt.ylabel('Energy (MHz)',fontsize=14)
        return

    def convert_evecs(self,basis,evecs=None,Normalize=True,verbose=True):
        if evecs is None:
            evecs = self.evecs0
        current_case = self.hunds_case
        new_case = basis
        if new_case == current_case:
            print('Eigenvectors already in {} basis'.format(current_case))
            return evecs
        inputt = self.q_numbers
        output = self.alt_q_numbers[new_case]
        if ('a' in new_case and 'bBJ' in current_case) or ('bBJ' in new_case and 'a' in current_case):
            basis_matrix = self.library.basis_changers['a_bBJ'](inputt,output)
        elif ('decoupled' in new_case and 'b' in current_case):
            basis_matrix = self.library.basis_changers['b_decoupled'](inputt,output)
        elif ('decoupled' in new_case and 'a' in current_case and 'J' not in new_case):
            intermediate = self.alt_q_numbers['bBJ']
            basis_matrix = self.library.basis_changers['b_decoupled'](intermediate,output)@self.library.basis_changers['a_bBJ'](inputt,intermediate)
        elif ('recoupled' in new_case and 'a' in current_case):
            intermediate1 = self.alt_q_numbers['bBJ']
            intermediate2 = self.alt_q_numbers['decoupled']
            basis_matrix = self.library.basis_changers['recouple_J'](intermediate2,output)@self.library.basis_changers['b_decoupled'](intermediate1,intermediate2)@self.library.basis_changers['a_bBJ'](inputt,intermediate1)
        elif ('a' in new_case and 'bBS' in current_case) or ('bBS' in new_case and 'a' in current_case):
            intermediate = self.alt_q_numbers['bBJ']
            if 'bBS' in current_case:
                basis_matrix = self.library.basis_changers['a_bBJ'](intermediate,output)@self.library.basis_changers['bBS_bBJ'](inputt,intermediate)
            else:
                basis_matrix = self.library.basis_changers['bBS_bBJ'](inputt,intermediate)@self.library.basis_changers['a_bBJ'](intermediate,output)
        converted_evecs = []
        for i in range(len(evecs)):
            converted_evecs.append(basis_matrix@evecs[i])
        converted_evecs = np.array(converted_evecs)
        if Normalize:
            for i,evec in enumerate(converted_evecs):
                converted_evecs[i]/=evec@evec
        if verbose:
            print('Successfully converted eigenvectors from {} to {}'.format(current_case,new_case))
        return converted_evecs

    def gen_state_str(self,vector_idx,evecs=None,basis=None,label_q=None,parity=False,single=False,thresh=0.01,show_coeff=True,new_line=False,round=None,frac=''):
        q_numbers = self.q_numbers
        if label_q == None:
            label_q = self.q_str
        if round is None:
            round=self.round
        if evecs is None:
            evecs = self.evecs0
        if basis is not None:
            if basis in self.hunds_case:
                pass
            elif 'decouple' in basis:
                q_numbers = self.alt_q_numbers['decoupled']
                if label_q == None:
                    label_q = list(self.alt_q_numbers['decoupled'])
                evecs = self.convert_evecs('decoupled',evecs=evecs,verbose=False)
            elif 'a' in basis:
                q_numbers = self.alt_q_numbers['aBJ']
                if label_q == None:
                    label_q = list(self.alt_q_numbers['aBJ'])
                evecs = self.convert_evecs('aBJ',evecs=evecs,verbose=False)
            elif 'bBJ' in basis:
                q_numbers = self.alt_q_numbers['bBJ']
                if label_q == None:
                    label_q = list(self.alt_q_numbers['bBJ'])
                evecs = self.convert_evecs('bBJ',evecs=evecs,verbose=False)
            elif 'recouple' in basis:
                q_numbers = self.alt_q_numbers['recouple_J']
                if label_q == None:
                    label_q = list(self.alt_q_numbers['recouple_J'])
            elif 'bBS' in basis:
                q_numbers = self.alt_q_numbers['bBS']
                if label_q == None:
                    label_q = list(self.alt_q_numbers['bBS'])
                evecs = self.convert_evecs('bBS',evecs=evecs,verbose=False)
        full_label = r''
        vector = deepcopy(evecs[vector_idx])
        vector[abs(vector)<thresh] = 0
        nonzero_idx = np.nonzero(vector)[0]
        if single:
            max_idx = abs(vector).argmax()
            nonzero_idx = [max_idx]
        i=0
        for i,index in enumerate(nonzero_idx):
            coeff = np.round(vector[index],round)
            sign = {True: '+', False: '-'}[coeff > 0]
            sign0 = {True: ' ', False: sign}[sign=='+']
            if show_coeff:
                label_str = r'$\,{}\,{}|'.format({True: sign0, False: sign}[i==0],abs(coeff))
            else:
                if single:
                    label_str = r'$|'
                else:
                    label_str = r'$\,{}|'.format({True: sign0, False: sign}[i==0])
                    if i==0:
                        label_str = r'$'+label_str[3:]
            if new_line:
                label_str = r'$'+label_str
            val = {q:q_numbers[q][index] for q in label_q}
            if parity:
                label_str+= r'{},'.format({1:'+',-1:'-'}[self.parities[vector_idx]])
            for q in label_q:
                _q = {True: '\u039B', False: q}[q=='L']
                _q = {True: '\u03A3', False: q}[q=='Sigma']
                _q = {True: '\u03A9', False: q}[q=='Omega']
                if (abs(val[q]) % 1) !=0:
                    label_str+=r'{}=\{}frac{{{}}}{{{}}},'.format(_q,frac,*val[q].as_integer_ratio())
                else:
                    label_str+=r'{}={},'.format(_q,int(val[q]))
            label_str = label_str[:-1]+r'\rangle\,$'
            if new_line:
                label_str+=r'$'
            full_label+=label_str
        return full_label

    def select_q(self,q_dict,evecs=None,parity=None):
        if evecs == None:
            evecs = self.evecs0
        idx = []
        for i in range(len(evecs)):
            match = True
            for q_string in q_dict:
                if not isinstance(q_dict[q_string],list):
                    q_dict[q_string] = [q_dict[q_string]]
                if self.q_numbers[q_string][np.argmax(evecs[i]**2)] in q_dict[q_string]:
                    match*=True
                else:
                    match*=False
                if parity is not None:
                    if parity=='+' and self.parities[i]<0:
                        match*=False
                    if parity=='-' and self.parities[i]>0:
                        match*=False
            if match==True:
                idx.append(i)
        return np.array(idx)

    def EB_grid(self,Ez,Bz,E_or_B_first = 'E',interp=False,output=False,evecs=False,PTV=False,EDM_or_MQM='EDM'):
        N_evals = len(self.evals0)
        evec_dim = len(self.evecs0[0])
        N_Bz = len(Bz)
        N_Ez = len(Ez)
        if '174' in self.iso_state or '40' in self.iso_state:
            self.PTV_type = 'EDM'
            H_PTV = self.library.PTV_builders[self.iso_state](self.q_numbers)
        elif '173' in self.iso_state or '171' in self.iso_state:
            self.PTV_type = EDM_or_MQM
            H_PTV = self.library.PTV_builders[self.iso_state](self.q_numbers, EDM_or_MQM)
        evals_EB = np.zeros((N_evals, N_Ez, N_Bz))
        if evecs:
            evecs_EB = np.zeros((N_evals, N_Ez, N_Bz,evec_dim))
        if PTV:
            PTV_EB = np.zeros((N_evals, N_Ez, N_Bz))
        if E_or_B_first == 'E':
            for i in range(N_Ez):
                evals_B, evecs_B = self.ZeemanMap(Bz, Ez_val = Ez[i],output=True)
                for j in range(N_Bz):
                    evals_EB[:,i,j] = evals_B[j]
                    if evecs:
                        evecs_EB[:,i,j,:] = evecs_B[j]
                    if PTV:
                        PTV_EB[:,i,j] = np.round(np.diagonal(evecs_B[j]@H_PTV@evecs_B[j].T),self.round)
                        #I think it has to be evecs.T from trial and error
        else:
            for i in range(N_Bz):
                evals_E, evecs_E = self.StarkMap(Ez, Bz_val = Bz[i], output=True)
                for j in range(N_Ez):
                    evals_EB[:,j,i] = evals_E[j]
                    if evecs:
                        evecs_EB[:,j,i,:] = evecs_E[j]
                    if PTV:
                        PTV_EB[:,j,i] = np.round(np.diagonal(evecs_E[j]@H_PTV@evecs_E[j].T),self.round)
        self.evals_EB = evals_EB
        if evecs:
            self.evecs_EB = evecs_EB
        if PTV:
            self.PTV_EB = PTV_EB
        if output and not (evecs or PTV):
            result = evals_EB
            return result
        elif output and (evecs or PTV):
            result = [evals_EB]
            if evecs:
                result.append(evecs_EB)
            if PTV:
                result.append(PTV_EB)
            return result
        else:
            return


    def generate_parities(self,evecs=None,ret=True):
        if evecs is None:
            evals,evecs = self.eigensystem(0,1e-6,set_attr=False)
        P_matrix = self.Parity_mat
        parities = []
        for evec in evecs:
            parity = np.round(evec@P_matrix@evec,0)
            parities.append(parity)
        self.parities = parities
        if ret:
            return parities
        return

def branching_ratios(Ground, Excited,Ez, Bz, Normalize=False):
    G_evals,G_evecs = Ground.eigensystem(Ez,Bz)
    G_qn = Ground.q_numbers
    E_evals,E_evecs = Excited.eigensystem(Ez,Bz)
    E_qn = Excited.q_numbers
    if 'a' not in Ground.hunds_case:
        G_evecs = Ground.convert_evecs('aBJ',Normalize=Normalize)
        G_qn = Ground.alt_q_numbers['aBJ']
    if 'a' not in Excited.hunds_case:
        E_evecs = Excited.convert_evecs('aBJ',Normalize=Normalize)
        E_qn = Excited.alt_q_numbers['aBJ']
    TDM_matrix = Excited.library.TDM_builders[Excited.iso_state](E_qn,G_qn)
    BR_matrix = (G_evecs@TDM_matrix@E_evecs.T)**2
    return BR_matrix

def XA_branching_ratios(X,A,Ez,Bz,Normalize=False): # must be in case a
    A.eigensystem(Ez,Bz)
    X.eigensystem(Ez,Bz)
    X_evecs_a = X.convert_evecs('aBJ',Normalize=Normalize)
    TDM_matrix = A.library.TDM_builders[A.iso_state](A.q_numbers,X.alt_q_numbers['aBJ'])
    BR_matrix = (X_evecs_a@TDM_matrix@A.evecs0.T)**2
    return BR_matrix

def Calculate_TDMs(p,Ground, Excited, Ez, Bz, Normalize=False):
    G_evals,G_evecs = Ground.eigensystem(Ez,Bz)
    G_qn = Ground.q_numbers
    E_evals,E_evecs = Excited.eigensystem(Ez,Bz)
    E_qn = Excited.q_numbers
    if 'a' not in Ground.hunds_case:
        G_evecs = Ground.convert_evecs('aBJ',Normalize=Normalize)
        G_qn = Ground.alt_q_numbers['aBJ']
    if 'a' not in Excited.hunds_case:
        E_evecs = Excited.convert_evecs('aBJ',Normalize=Normalize)
        E_qn = Excited.alt_q_numbers['aBJ']
    TDM_matrix = Excited.library.TDM_p_builders[Excited.iso_state](p,E_qn,G_qn)
    TDM_p = (G_evecs@TDM_matrix@E_evecs.T)
    return TDM_p


def state_ordering(evecs_old,evecs_new,step=None,round=8):
    if step is None:
        step = 1
    overlap = abs(evecs_old@evecs_new.T)     #Essentially a matrix of the fidelities: |<phi|psi>|
    #calculate trace distance
    # for o in overlap:
    #     for _o in o:
    #         if (_o>1):
    #             print('OVERLAP BIGGER THAN 1', _o)
    dist = abs(1-overlap)/step**2
    # ordering = np.array([dist[i,:].argmin() for i in range(len(evecs_old))])  #python
    ordering = np.argmin(dist,axis=1) #numpy
    return ordering

def order_eig(evals,evecs):
    order = np.argsort(evals)
    evecs_ordered =evecs[order,:]
    evals_ordered = evals[order]
    return evals_ordered,evecs_ordered

def diagonalize_batch(matrix_array,method='torch',round=10):
    if method == 'numpy':
        w,v = npLA.eigh(matrix_array)
    elif method == 'scipy':
        w,v = sciLA.eigh(matrix_array)
    elif method == 'torch':
        w,v = torch.linalg.eigh(torch.from_numpy(matrix_array))
        w = w.numpy()
        v = v.numpy()
    evals_batch = w
    evecs_batch= np.transpose(v,[0,2,1])
    return evals_batch,evecs_batch


def diagonalize(matrix,method='torch',order=False, Normalize=False,round=10):
    if method == 'numpy':
        w,v = npLA.eigh(matrix)
    elif method == 'scipy':
        w,v = sciLA.eigh(matrix)
    elif method == 'torch':
        w,v = torch.linalg.eigh(torch.from_numpy(matrix))
        w = w.numpy()
        v = v.numpy()
    evals = w
    evecs = v.T
    # if Normalize:
    #     for i,evec in enumerate(evecs):
    #         evecs[i]/=evec@evec
    # if order:
    #     idx_order = np.argsort(evals)
    #     evecs =evecs[idx_order,:]
    #     evals = evals[idx_order]
    return evals,evecs
