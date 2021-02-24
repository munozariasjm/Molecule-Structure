import YbOH_matrix_elements as me
import YbOH_quantum_numbers as qn
import YbOH_hamiltonian_builders as ham
from YbOH_parameters import all_params
from functools import partial

class YbOH_Library(object):

    '''
    This class is primarily used as a library of functions and properties used
    to model the structure of YbOH. The class attributes are:

    matrix_elements: dictionary of relevant matrix elements, indexed by the
    YbOH state and isotope (ie, '174X000')

    cases: dictionary of the Hund's cases used for each state. The matrix
    elements are calculated in this basis.

    H_builders: dictionary of functions that construct the Hamiltonian matrix
    for a given isotope and state.

    PTV_builders: dictionary of functions that construct Hamiltonians for
    PT-violating effects.

    q_number_builders: dictionary of functions that construct the quantum
    numbers for a given state and isotope.
    The quantum numbers are represented as a dictionary, with the dict key
    giving the angular momentum name, and the dict value containing an array of
    the eigenvalues of each basis vector.
    For example, {'J': [0,1,1,1], 'M': [0,-1,0,1]}

    parameters: dictionary of the Hamiltonian constants associated with each
    state and isotope.

    Lambda: projection of angular momentum on the internuclear axis, excluding
    electron spin.
    '''

    def __init__(self,I_spins):
        self.parameters = all_params
        self.matrix_elements = self.collect_all_matrix_elements(I_spins)
        self.cases = self.collect_all_cases()
        self.H_builders = self.collect_all_H_builders()
        self.PTV_builders = self.collect_all_PTV_builders()
        self.q_number_builders = self.collect_all_q_number_builders()
        self.Lambda = self.collect_all_Lambda()
        self.basis_changers = self.collect_change_basis()
        self.TDM_builders = self.collect_TDM()
        self.alt_q_number_builders = self.collect_alt_q()


    def collect_all_cases(self):
        all_cases={
        '174X000': 'bBJ',
        '174X010': 'bBJ',
        '173X000': 'bBS',
        '173X010': 'bBS',
        '174A000': 'aBJ',
        '173A000': 'aBJ_seq'
        }
        return all_cases

    def collect_all_Lambda(self):
        all_Lambda={
        '174X000': 0,
        '174X010': 1,
        '173X000': 0,
        '173X010': 1,
        '174A000': 1,
        '173A000': 1
        }
        return all_Lambda


    def collect_all_matrix_elements(self,I_spins):
        iH = I_spins[-1]
        IYb = I_spins[0]

        bBJ_174X_matrix_elements={
        # Fine Structure
        'N^2': me.Rot_bBJ,                 # N^2 Rotation
        'N.S': me.SR_bBJ,                  # N.S Spin Rotation
        'l-doubling': me.lD_bBJ,           # Effective l doubling

        # Hydrogen Hyperfine
        'I.S': me.IS_bBJ,                  # I.S Fermi Contact Interaction
        'T2_0(I,S)': me.T2IS_bBJ,          # I S dipolar interaction
        'Iz': me.Iz_bBJ,                   # I.n projection of I on internuclear axis n
        'Sz': me.Sz_bBJ,                   # S.n projection of S on internuclear axis n

        # External Fields
        'ZeemanZ': me.ZeemanZ_bBJ,         # Zeeman interaction with lab z magnetic field
        'StarkZ': me.StarkZ_bBJ            # Stark interaction with lab z electric field
        }
        for term,element in bBJ_174X_matrix_elements.items():       #iterate through, substitute hydrogen proton value
            bBJ_174X_matrix_elements[term] = partial(element,I=iH)




        bBS_173X_matrix_elements={
        # Fine Structure
        'N^2': me.Rot_bBS,                 # N^2 Rotation
        'N.S': me.SR_bBS,                  # N.S Spin Rotation
        'l-doubling': me.lD_bBS,           # Effective l doubling

        # Yb Hyperfine
        'IYb.S': me.ISYb_bBS,              # IYb.S Fermi Contact Interaction
        'T2_0(IYb,S)': me.T2ISYb_bBS,      # q=0 molecule frame electron-IYb dipole-dipole interaction (sqrt6 * (3IzSz - I.S))
        'T2_0(IYb^2)': me.T2QYb_bBS,       # q=0 molecule frame Yb electric quadrupole interaction (sqrt6 * (3Iz^2 - I^2))

        # Hydrogen Hyperfine
        'IH.S': me.ISH_bBS,                # IH.S Fermi Contact Interaction for Hydrogen
        'T2_0(IH,S)': me.T2ISH_bBS,        # q=0 molecule frame electron-IH dipole-dipole interaction (sqrt6 * (3IzSz - I.S))

        # External Fields
        'ZeemanZ': me.ZeemanZ_bBS,         # Zeeman interaction with lab z magnetic field
        'StarkZ': me.StarkZ_bBS            # Stark interaction with lab z electric field
        }
        for term, element in bBS_173X_matrix_elements.items():      #substitue nuclear spin values
             bBS_173X_matrix_elements[term] = partial(element,iH=iH,I=IYb)

        aBJ_174A_matrix_elements={
        # Fine Structure
        'N^2': me.Rot_174_aBJ,                 # N^2 Rotation
        'SO': me.SO_174_aBJ,                   # Spin-Orbit
        'Lambda-Doubling': me.LambdaDoubling_174_aBJ,          # Effective Lambda doubling

        # Hydrogen Hyperfine
        'I.S': me.IS_174_aBJ,                   # I.S Fermi Contact Interaction
        'IzSz': me.IzSz_174_aBJ,                # I.n*S.n projection of I and S on internuclear axis n
        'T2_0(IS)': me.T2IS_174_aBJ,            # q=0 molecule frame electron-IH dipole-dipole interaction (sqrt6 * (3IzSz - I.S))

        # External Fields
        'ZeemanLZ': me.ZeemanLZ_174_aBJ,
        'ZeemanSZ': me.ZeemanSZ_174_aBJ,
        'ZeemanParityZ': me.ZeemanParityZ_174_aBJ,
        'StarkZ': me.StarkZ_174_aBJ,            # Stark interaction with lab z electric field
        }
        for term, element in aBJ_174A_matrix_elements.items():      #substitue nuclear spin values
             aBJ_174A_matrix_elements[term] = partial(element,I=iH)


        aBJ_173A_matrix_elements={
        # Fine Structure
        'N^2': me.Rot_173_aBJ,                 # N^2 Rotation
        'Lamdba-Doubling': me.LambdaDoubling_173_aBJ,

        # Yb Hypeerfine
        'IzLz_Yb': me.ILYb_173_aBJ,
        'T2_2(IS)_Yb': me.T2q2_ISYb_173_aBJ,
        'T2_0(IS)_Yb': me.T2q0_IYb_173_aBJ,


        # Hydrogen Hyperfine
        # 'IS_H': ,                   # I.S Fermi Contact Interaction
        # 'Iz_H': ,                   # I.n projection of I on internuclear axis n
        # 'Sz_H': ,                   # S.n projection of S on internuclear axis n
        #
        # # External Fields
        # 'ZeemanZ': me.ZeemanZ_aBJ,         # Zeeman interaction with lab z magnetic field
        # 'StarkZ': me.StarkZ_a            # Stark interaction with lab z electric field
        }
        for term, element in aBJ_173A_matrix_elements.items():      #substitue nuclear spin values
             aBJ_173A_matrix_elements[term] = partial(element,iH=iH,I=IYb)

        all_matrix_elements={
        '174X000': bBJ_174X_matrix_elements,
        '174X010': bBJ_174X_matrix_elements,
        '173X000': bBS_173X_matrix_elements,
        '173X010': bBS_173X_matrix_elements,
        '174A000': aBJ_174A_matrix_elements,
        '173A000': aBJ_173A_matrix_elements
        }

        return all_matrix_elements

    def collect_all_H_builders(self):
        H_builders = {
            '174X000': ham.H_174X,
            '174X010': ham.H_174X,
            '173X000': ham.H_173X,
            '173X010': ham.H_173X,
            '174A000': ham.H_174A,
            '173A000': ham.H_173A
            }
        for key in H_builders:
            H_builders[key] = partial(H_builders[key],params = self.parameters[key],matrix_elements = self.matrix_elements[key])
        return H_builders

    def collect_all_q_number_builders(self):
        q_number_builders = {
            '174X000': partial(qn.q_numbers_bBJ, Lambda=0),
            '174X010': partial(qn.q_numbers_bBJ, Lambda=1),
            '173X000': partial(qn.q_numbers_bBS, Lambda=0),
            '173X010': partial(qn.q_numbers_bBS, Lambda=1),
            '174A000': partial(qn.q_numbers_aBJ, Lambda=1,Omega_values=[1/2]),
            '173A000': partial(qn.q_numbers_aBJ, Lambda=1,Omega_values=[1/2]),
            '174aBJ': qn.q_numbers_aBJ,
            '174bBJ': qn.q_numbers_bBJ,
        }
        return q_number_builders

    def collect_all_PTV_builders(self):
        PTV_builders = {
            '174X000': ham.build_PTV_bBJ,
            '174X010': ham.build_PTV_bBJ,
            '173X000': ham.build_PTV_bBS,
            '173X010': ham.build_PTV_bBS
        }
        return PTV_builders

    def collect_change_basis(self):
        all_change_basis = {
            'a_b': ham.convert_ab,
            'b_decoupled': ham.decouple_b
        }
        return all_change_basis

    def collect_TDM(self):
        all_TDM = {
            '174A000': partial(ham.build_TDM_174_aBJ,TDM_matrix_element=me.TransitionDipole_174_aBJ)
        }
        return all_TDM

    def collect_alt_q(self):
        alt_q_builders = {
            '174X000': {
                'aBJ': partial(qn.q_numbers_aBJ, Lambda=0),
                'decoupled': partial(qn.q_numbers_decoupled, Lambda=0)
                },
            '174X010': {
                'aBJ': partial(qn.q_numbers_aBJ, Lambda=1,Omega_values=[1/2,3/2]),
                'decoupled': partial(qn.q_numbers_decoupled, Lambda=1)
                },
            '173X000': {
                'aBJ': partial(qn.q_numbers_aBJ, Lambda=0),
                'decoupled': partial(qn.q_numbers_decoupled,Lambda=0)
            },
            '173X010': {
                'aBJ': partial(qn.q_numbers_aBJ, Lambda=1,Omega_values=[1/2,3/2]),
                'decoupled': partial(qn.q_numbers_decoupled, Lambda=1)
            },
            '174A000': {
                'bBS': partial(qn.q_numbers_bBJ, Lambda=1),
                'decoupled': partial(qn.q_numbers_decoupled, Lambda=1)
            },
            '173A000': {
                'bBS': partial(qn.q_numbers_bBS, Lambda=1),
                'decoupled': partial(qn.q_numbers_decoupled, Lambda=1)
            }
        }
        return alt_q_builders
