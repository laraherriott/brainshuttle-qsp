from .parameters import NoAbParameters


class NoAbModel:
    """Class representing the ODE model with no antibody and two Abeta
    populations. These can be used to represent Abeta42 and Abeta40, or
    leucine labeled and unlabeled Abeta.

    """
    def __init__(self, leucine=None, x=None, healthy=False, time_points=None):
        """Constructor Method.

        Parameters
        ----------
        leucine : List
            Nested list, [List, List] representing the proportion of leucine
            labeled at each time point in the plasma and CNS respectively
        x : List
            Input parameters to override those provided in
            class:NoAbParameters if desired
        healthy : Boolean
            Whether the base model represents a healthy individual
            or an AD partient at the MCI stage
        time_points : Int
            The number of time points overwhich the simulation will be run

        """
        if time_points is None:
            raise ValueError("The total number of timepoints for\
                              the simulation must be provided")
        if x is None:
            self.params = NoAbParameters()
        elif healthy:
            self.general = NoAbParameters()
            self.general.k_in = 10**x[0]
            self.general.k_peripheral_production = 10**x[1]
            self.general.k_olig_inc = 10**x[2]
            self.general.k_olig_sep = 10**x[3]
            self.general.k_olig_inc_ext = 10**x[4]
            self.general.k_olig_sep_ext = 10**x[5]
            self.general.k_clear_Abeta_plasma = 10**x[6]
            self.general.k_clear_oligomer_plasma = 10**x[7]
            self.general.k_clear_Abeta_brain = 10**x[8]
            self.general.k_clear_oligomer_brain = 10**x[9]
            self.general.k_monomer_plasma_brain = 10**x[10]
            self.general.k_oligomer_plasma_brain = 10**x[11]
            self.general.k_monomer_brain_plasma = 10**x[12]
            self.general.k_oligomer_brain_plasma = 10**x[13]
            self.general.k_monomer_brain_csf = 10**x[14]
            self.general.k_oligomer_brain_csf = 10**x[15]
            self.general.k_monomer_csf_brain = 10**x[16]
            self.general.k_oligomer_csf_brain = 10**x[17]
            self.general.k_monomer_csf_plasma = 10**x[18]
            self.general.k_oligomer_csf_plasma = 10**x[19]
            self.general.k_plaque_inc = 0
            self.general.k_plaque_sep = 0
            self.general.k_clear_plaque = 0
            self.params = self.general
        else:
            self.general = NoAbParameters()
            self.general.k_in = 10**x[0]
            self.general.k_peripheral_production = 10**x[1]
            self.general.k_olig_inc = 10**x[2]
            self.general.k_olig_sep = 10**x[3]
            self.general.k_olig_inc_ext = 10**x[4]
            self.general.k_olig_sep_ext = 10**x[5]
            self.general.k_plaque_inc = 10**x[6]
            self.general.k_plaque_sep = 10**x[7]
            self.general.k_clear_Abeta_plasma = 10**x[8]
            self.general.k_clear_oligomer_plasma = 10**x[9]
            self.general.k_clear_Abeta_brain = 10**x[10]
            self.general.k_clear_oligomer_brain = 10**x[11]
            self.general.k_monomer_plasma_brain = 10**x[12]
            self.general.k_oligomer_plasma_brain = 10**x[13]
            self.general.k_monomer_brain_plasma = 10**x[14]
            self.general.k_oligomer_brain_plasma = 10**x[15]
            self.general.k_monomer_brain_csf = 10**x[16]
            self.general.k_oligomer_brain_csf = 10**x[17]
            self.general.k_monomer_csf_brain = 10**x[18]
            self.general.k_oligomer_csf_brain = 10**x[19]
            self.general.k_monomer_csf_plasma = 10**x[20]
            self.general.k_oligomer_csf_plasma = 10**x[21]
            self.general.k_clear_plaque = 10**x[22]
            self.params = self.general

        if leucine is not None:
            self.f_leu = leucine
        else:
            self.f_leu = [[0]*time_points, [0]*time_points]

        if healthy:
            self.type = 'healthy'
        else:
            self.type = 'no_ab'

    def equations(self, t, y):
        """Ordinary Differential Equations used by scipy.integrate.solve_ivp()

        Parameters
        ----------
        t : float
            Current time
        y : List
            Current values for each variable

        Returns
        ----------
        dYdt : List
            Change in each variable at current time

        """
        unlabeled_production = (1-self.f_leu[1][int(t)])*self.params.k_in
        unlabeled_peripheral_production = (1-self.f_leu[0][int(t)])*self.\
            params.k_peripheral_production

        labeled_production = (self.f_leu[1][int(t)])*self.params.k_in
        labeled_peripheral_production = (self.f_leu[0][int(t)])*self.\
            params.k_peripheral_production

        Vplasma = 3.5
        Vcsf = 0.15
        Visf = 0.3

        brain_monomer42 = y[0]
        brain_monomer40 = y[1]

        brain_oligomer42 = y[2]
        brain_oligomer40 = y[3]

        brain_plaque42 = y[4]
        brain_plaque40 = y[5]

        plasma_monomer42 = y[6]
        plasma_monomer40 = y[7]

        plasma_oligomer42 = y[8]
        plasma_oligomer40 = y[9]

        csf_monomer42 = y[10]
        csf_monomer40 = y[11]

        csf_oligomer42 = y[12]
        csf_oligomer40 = y[13]

        d_brain_monomer42 = (unlabeled_production
                             - self.params.k_olig_inc*brain_monomer42
                             + self.params.k_olig_sep*brain_oligomer42
                             - self.params.k_clear_Abeta_brain*brain_monomer42
                             - self.params.k_monomer_brain_plasma*brain_monomer42
                             + self.params.k_monomer_plasma_brain*(Vplasma/Visf)*plasma_monomer42
                             - self.params.k_monomer_brain_csf*brain_monomer42
                             + self.params.k_monomer_csf_brain*(Vcsf/Visf)*csf_monomer42)

        d_brain_monomer40 = (labeled_production
                             - self.params.k_olig_inc*brain_monomer40
                             + self.params.k_olig_sep*brain_oligomer40
                             - self.params.k_clear_Abeta_brain*brain_monomer40
                             - self.params.k_monomer_brain_plasma*brain_monomer40
                             + self.params.k_monomer_plasma_brain*(Vplasma/Visf)*plasma_monomer40
                             - self.params.k_monomer_brain_csf*brain_monomer40
                             + self.params.k_monomer_csf_brain*(Vcsf/Visf)*csf_monomer40)

        d_brain_oligomer40 = (self.params.k_olig_inc*brain_monomer40
                              - self.params.k_olig_sep*brain_oligomer40
                              - self.params.k_plaque_inc*brain_oligomer40
                              + self.params.k_plaque_sep*brain_plaque40
                              - self.params.k_clear_oligomer_brain*brain_oligomer40
                              - self.params.k_oligomer_brain_plasma*brain_oligomer40
                              + self.params.k_oligomer_plasma_brain*(Vplasma/Visf)*plasma_oligomer40
                              - self.params.k_oligomer_brain_csf*brain_oligomer40
                              + self.params.k_oligomer_csf_brain*(Vcsf/Visf)*csf_oligomer40)

        d_brain_oligomer42 = (self.params.k_olig_inc*brain_monomer42
                              - self.params.k_olig_sep*brain_oligomer42
                              - self.params.k_plaque_inc*brain_oligomer42
                              + self.params.k_plaque_sep*brain_plaque42
                              - self.params.k_clear_oligomer_brain*brain_oligomer42
                              - self.params.k_oligomer_brain_plasma*brain_oligomer42
                              + self.params.k_oligomer_plasma_brain*(Vplasma/Visf)*plasma_oligomer42
                              - self.params.k_oligomer_brain_csf*brain_oligomer42
                              + self.params.k_oligomer_csf_brain*(Vcsf/Visf)*csf_oligomer42)

        d_brain_plaque40 = (self.params.k_plaque_inc*brain_oligomer40
                            - self.params.k_plaque_sep*brain_plaque40
                            - self.params.k_clear_plaque*brain_plaque40)
        d_brain_plaque42 = (self.params.k_plaque_inc*brain_oligomer42
                            - self.params.k_plaque_sep*brain_plaque42
                            - self.params.k_clear_plaque*brain_plaque42)

        d_plasma_monomer42 = (unlabeled_peripheral_production
                              + self.params.k_monomer_brain_plasma*(Visf/Vplasma)*brain_monomer42
                              - self.params.k_monomer_plasma_brain*plasma_monomer42
                              + self.params.k_monomer_csf_plasma*(Vcsf/Vplasma)*csf_monomer42
                              - self.params.k_monomer_plasma_csf*plasma_monomer42
                              - self.params.k_olig_inc_ext*plasma_monomer42
                              + self.params.k_olig_sep_ext*plasma_oligomer42
                              - self.params.k_clear_Abeta_plasma*plasma_monomer42)

        d_plasma_monomer40 = (labeled_peripheral_production
                              + self.params.k_monomer_brain_plasma*(Visf/Vplasma)*brain_monomer40
                              - self.params.k_monomer_plasma_brain*plasma_monomer40
                              + self.params.k_monomer_csf_plasma*(Vcsf/Vplasma)*csf_monomer40
                              - self.params.k_monomer_plasma_csf*plasma_monomer40
                              - self.params.k_olig_inc_ext*plasma_monomer40
                              + self.params.k_olig_sep_ext*plasma_oligomer40
                              - self.params.k_clear_Abeta_plasma*plasma_monomer40)

        d_plasma_oligomer40 = (self.params.k_oligomer_brain_plasma*(Visf/Vplasma)*brain_oligomer40
                               - self.params.k_oligomer_plasma_brain*plasma_oligomer40
                               + self.params.k_oligomer_csf_plasma*(Vcsf/Vplasma)*csf_oligomer40
                               - self.params.k_oligomer_plasma_csf*plasma_oligomer40
                               - self.params.k_clear_oligomer_plasma*plasma_oligomer40
                               + self.params.k_olig_inc_ext*plasma_monomer40
                               - self.params.k_olig_sep_ext*plasma_oligomer40)

        d_plasma_oligomer42 = (self.params.k_oligomer_brain_plasma*(Visf/Vplasma)*brain_oligomer42
                               - self.params.k_oligomer_plasma_brain*plasma_oligomer42
                               + self.params.k_oligomer_csf_plasma*(Vcsf/Vplasma)*csf_oligomer42
                               - self.params.k_oligomer_plasma_csf*plasma_oligomer42
                               - self.params.k_clear_oligomer_plasma*plasma_oligomer42
                               + self.params.k_olig_inc_ext*plasma_monomer42
                               - self.params.k_olig_sep_ext*plasma_oligomer42)

        d_csf_monomer42 = (self.params.k_monomer_brain_csf*(Visf/Vcsf)*brain_monomer42
                           - self.params.k_monomer_csf_brain*csf_monomer42
                           - self.params.k_monomer_csf_plasma*csf_monomer42
                           + self.params.k_monomer_plasma_csf*(Vplasma/Vcsf)*plasma_monomer42
                           - self.params.k_olig_inc_ext*csf_monomer42
                           + self.params.k_olig_sep_ext*csf_oligomer42
                           - self.params.k_clear_Abeta_csf*csf_monomer42)

        d_csf_monomer40 = (self.params.k_monomer_brain_csf*(Visf/Vcsf)*brain_monomer40
                           - self.params.k_monomer_csf_brain*csf_monomer40
                           - self.params.k_monomer_csf_plasma*csf_monomer40
                           + self.params.k_monomer_plasma_csf*(Vplasma/Vcsf)*plasma_monomer40
                           - self.params.k_olig_inc_ext*csf_monomer40
                           + self.params.k_olig_sep_ext*csf_oligomer40
                           - self.params.k_clear_Abeta_csf*csf_monomer40)

        d_csf_oligomer40 = (self.params.k_oligomer_brain_csf*(Visf/Vcsf)*brain_oligomer40
                            - self.params.k_oligomer_csf_brain*csf_oligomer40
                            - self.params.k_oligomer_csf_plasma*csf_oligomer40
                            + self.params.k_oligomer_plasma_csf*(Vplasma/Vcsf)*plasma_oligomer40
                            - self.params.k_clear_oligomer_csf*csf_oligomer40
                            + self.params.k_olig_inc_ext*csf_monomer40
                            - self.params.k_olig_sep_ext*csf_oligomer40)

        d_csf_oligomer42 = (self.params.k_oligomer_brain_csf*(Visf/Vcsf)*brain_oligomer42
                            - self.params.k_oligomer_csf_brain*csf_oligomer42
                            - self.params.k_oligomer_csf_plasma*csf_oligomer42
                            + self.params.k_oligomer_plasma_csf*(Vplasma/Vcsf)*plasma_oligomer42
                            - self.params.k_clear_oligomer_csf*csf_oligomer42
                            + self.params.k_olig_inc_ext*csf_monomer42
                            - self.params.k_olig_sep_ext*csf_oligomer42)

        dYdt = [d_brain_monomer42, d_brain_monomer40,
                d_brain_oligomer42, d_brain_oligomer40,
                d_brain_plaque42, d_brain_plaque40,
                d_plasma_monomer42, d_plasma_monomer40,
                d_plasma_oligomer42, d_plasma_oligomer40,
                d_csf_monomer42, d_csf_monomer40,
                d_csf_oligomer42, d_csf_oligomer40]

        return dYdt
