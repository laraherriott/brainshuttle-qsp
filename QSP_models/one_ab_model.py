import math

from .parameters import Lecanemab, NoAbParameters


class OneAbModel:
    """Class representing the ODE model with one antibody.

    """
    def __init__(self, dose, bodyweight, gpm, final_dose=None, dose_interval=None,
                 dose_times=None,
                 x=None, y=None, z=None, healthy=False):
        """Constructor Method.

        Parameters
        ----------
        dose : Numeric or List
            Dose in mg/kg. If list represents size of each dose
            over whole simulation. If numeric same dose will be
            given at each time point.
        bodyweight : Numeric
            Representative bodyweight of the population in kg
        gpm : Numeric
            Grams per mole of antibody
        final dose : Int
            Time of final dose in dosing schedule
        dose_interval : Int
            Time between repeated doses in dosing schedule
        dose_times : List
            List of dose times in months if not regular schedule
        x : List
            Input parameters to override those provided in
            class:NoAbParameters if desired
        y : List
            Input parameters to override those provided in
            class:Lecanemab if desired
        z : List
            Input parameters to override default KD values for
            monomer, oligomer, and plaque, if desired
        healthy : Boolean
            Whether the base model represents a healthy individual
            or an AD partient at the MCI stage or a cognitively healthy individual

        """
        if x is None:
            self.general = NoAbParameters()
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

        if y is None:
            self.params = Lecanemab()
        elif healthy:
            self.params = Lecanemab()
            self.params.k_mAb_plasma_brain = 10**y[0]
            self.params.k_mAb_brain_csf = 10**y[1]
            self.params.k_mAb_csf_brain = 10**y[2]
            self.params.k_mAb_csf_plasma = 10**y[3]
            self.params.k_mAb_plasma_csf = 10**y[4]
            self.params.plasma_clearance = 10**y[5]
            self.params.k_mAb_brain_plasma = 10**y[6]
            self.params.brain_clearance = 0
            self.params.k_ADCP = 0
            self.params.scale = 0
        else:
            self.params = Lecanemab()
            self.params.k_mAb_plasma_brain = 10**y[0]
            self.params.k_mAb_brain_csf = 10**y[1]
            self.params.k_mAb_csf_brain = 10**y[2]
            self.params.k_mAb_csf_plasma = 10**y[3]
            self.params.k_mAb_plasma_csf = 10**y[4]
            self.params.plasma_clearance = 10**y[5]
            self.params.k_mAb_brain_plasma = 10**y[6]
            self.params.brain_clearance = 0
            self.params.k_ADCP = 10**y[7]
            self.params.scale = y[8]

        if z is not None:
            self.params.offma0 = z[0]*0.001
            self.params.offma1 = z[1]*0.001
            self.params.offma2 = z[2]*0.001

        if healthy:
            self.type = 'one_ab_healthy'
        else:
            self.type = 'one_ab'

        if dose_times is None:
            if final_dose is None or dose_interval is None:
                raise ValueError('either an explicit dose list for the duration of the simulation, or a dose interval and final time must be provided')
            else:
                self.dose_list = []
                i = 0
                while i <= final_dose:
                    self.dose_list.append(i)
                    i += dose_interval
        else:
            if isinstance(dose_times, list):
                self.dose_list = dose_times
            else:
                raise TypeError('if providing explicity list of dose times, must be type list')

        if isinstance(dose, float) or isinstance(dose, int):
            self.dose = [dose]*len(self.dose_list)
        elif isinstance(dose, list):
            if len(dose) != len(self.dose_list):
                raise ValueError('dose list length dose not match the\
                                 number of doses')
            else:
                self.dose = dose
        else:
            raise TypeError('dose must be int, float or list')

        self.bw = bodyweight
        self.gpm = gpm

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
        Vplasma = 3.5
        Vcsf = 0.15
        Visf = 0.3

        brain_ab = y[0]
        brain_monomer = y[1]
        brain_oligomer = y[2]
        brain_plaque = y[3]

        plasma_ab = y[4]
        plasma_monomer = y[5]
        plasma_oligomer = y[6]

        csf_ab = y[7]
        csf_monomer = y[8]
        csf_oligomer = y[9]

        brain_monomer_mAb = y[10]
        brain_oligomer_mAb = y[11]
        brain_plaque_mAb = y[12]

        plasma_monomer_mAb = y[13]
        plasma_oligomer_mAb = y[14]
        csf_monomer_mAb = y[15]
        csf_oligomer_mAb = y[16]

        d_plasma_mAb = (self.dosefn(self.gpm, self.dose_list, t)
                        - self.params.plasma_clearance*plasma_ab
                        + self.params.k_mAb_csf_plasma*(Vcsf/Vplasma)*csf_ab
                        - self.params.k_mAb_plasma_csf*plasma_ab
                        + self.params.k_mAb_brain_plasma*(Visf/Vplasma)*brain_ab
                        - self.params.k_mAb_plasma_brain*plasma_ab
                        - self.params.onPP*plasma_ab*plasma_monomer
                        + self.params.offma0*plasma_monomer_mAb
                        - self.params.onPP*plasma_ab*plasma_oligomer
                        + self.params.offma1*plasma_oligomer_mAb)

        d_brain_mAb = (- self.params.brain_clearance*brain_ab
                       - self.params.k_mAb_brain_plasma*brain_ab
                       - self.params.k_mAb_brain_csf*brain_ab
                       + self.params.k_mAb_csf_brain*(Vcsf/Visf)*csf_ab
                       + self.params.k_mAb_plasma_brain*(Vplasma/Visf)*plasma_ab
                       - self.params.onPP*brain_ab*brain_monomer
                       + self.params.offma0*brain_monomer_mAb
                       - self.params.onPP*brain_ab*brain_oligomer
                       + self.params.offma1*brain_oligomer_mAb
                       - self.params.onPP*brain_ab*brain_plaque
                       + self.params.offma2*brain_plaque_mAb)

        d_csf_mAb = (- self.params.k_mAb_csf_plasma*csf_ab
                     + self.params.k_mAb_brain_csf*(Visf/Vcsf)*brain_ab
                     + self.params.k_mAb_plasma_csf*(Vplasma/Vcsf)*plasma_ab
                     - self.params.k_mAb_csf_brain*csf_ab
                     - self.params.onPP*csf_ab*csf_monomer
                     + self.params.offma0*csf_monomer_mAb
                     - self.params.onPP*csf_ab*csf_oligomer
                     + self.params.offma1*csf_oligomer_mAb)

        d_brain_monomer = (self.general.k_in
                           - self.general.k_olig_inc*brain_monomer
                           + self.general.k_olig_sep*brain_oligomer
                           - self.general.k_clear_Abeta_brain*brain_monomer
                           - self.general.k_monomer_brain_plasma*brain_monomer
                           + self.general.k_monomer_plasma_brain*(Vplasma/Visf)*plasma_monomer
                           - self.general.k_monomer_brain_csf*brain_monomer
                           + self.general.k_monomer_csf_brain*(Vcsf/Visf)*csf_monomer
                           - self.params.onPP*brain_ab*brain_monomer
                           + self.params.offma0*brain_monomer_mAb)

        d_brain_monomer_mAb = (self.params.onPP*brain_ab*brain_monomer
                               - self.params.offma0*brain_monomer_mAb
                               - self.params.brain_clearance*brain_monomer_mAb
                               - self.params.k_mAb_brain_plasma*brain_monomer_mAb
                               - self.params.k_mAb_brain_csf*brain_monomer_mAb
                               + self.params.k_mAb_csf_brain*(Vcsf/Visf)*csf_monomer_mAb
                               + self.params.k_mAb_plasma_brain*(Vplasma/Visf)*plasma_monomer_mAb)

        d_brain_oligomer = (self.general.k_olig_inc*brain_monomer
                            - self.general.k_olig_sep*brain_oligomer
                            - self.general.k_plaque_inc*brain_oligomer
                            + self.general.k_plaque_sep*brain_plaque
                            - self.general.k_clear_oligomer_brain*brain_oligomer
                            - self.general.k_oligomer_brain_plasma*brain_oligomer
                            + self.general.k_oligomer_plasma_brain*(Vplasma/Visf)*plasma_oligomer
                            - self.general.k_oligomer_brain_csf*brain_oligomer
                            + self.general.k_oligomer_csf_brain*(Vcsf/Visf)*csf_oligomer
                            - self.params.onPP*brain_ab*brain_oligomer
                            + self.params.offma1*brain_oligomer_mAb)

        d_brain_oligomer_mAb = (self.params.onPP*brain_ab*brain_oligomer
                                - self.params.offma1*brain_oligomer_mAb
                                - self.params.brain_clearance*brain_oligomer_mAb
                                - self.params.k_ADCP*brain_oligomer_mAb
                                - self.params.k_mAb_brain_plasma*brain_oligomer_mAb
                                - self.params.k_mAb_brain_csf*brain_oligomer_mAb
                                + self.params.k_mAb_csf_brain*(Vcsf/Visf)*csf_oligomer_mAb
                                + self.params.k_mAb_plasma_brain*(Vplasma/Visf)*plasma_oligomer_mAb)

        coclearance = self.params.k_ADCP*self.params.scale*brain_plaque_mAb*brain_plaque

        d_brain_plaque = (self.general.k_plaque_inc*brain_oligomer
                          - self.general.k_plaque_sep*brain_plaque
                          - self.params.onPP*brain_ab*brain_plaque
                          + self.params.offma2*brain_plaque_mAb
                          - self.general.k_clear_plaque*brain_plaque
                          - coclearance/(brain_plaque+0.001))

        d_brain_plaque_mAb = (self.params.onPP*brain_ab*brain_plaque
                              - self.params.offma2*brain_plaque_mAb
                              - (self.params.k_ADCP)*brain_plaque_mAb)

        d_plasma_monomer = (self.general.k_peripheral_production
                            + self.general.k_monomer_brain_plasma*(Visf/Vplasma)*brain_monomer
                            - self.general.k_monomer_plasma_brain*plasma_monomer
                            + self.general.k_monomer_csf_plasma*(Vcsf/Vplasma)*csf_monomer
                            - self.params.onPP*plasma_ab*plasma_monomer
                            + self.params.offma0*plasma_monomer_mAb
                            - self.general.k_clear_Abeta_plasma*plasma_monomer
                            - self.general.k_olig_inc_ext*plasma_monomer
                            + self.general.k_olig_sep_ext*plasma_oligomer)

        d_plasma_monomer_mAb = (self.params.onPP*plasma_ab*plasma_monomer
                                - self.params.offma0*plasma_monomer_mAb
                                - self.params.plasma_clearance*plasma_monomer_mAb
                                + self.params.k_mAb_csf_plasma*(Vcsf/Vplasma)*csf_monomer_mAb
                                + self.params.k_mAb_brain_plasma*(Visf/Vplasma)*brain_monomer_mAb
                                - self.params.k_mAb_plasma_csf*plasma_monomer_mAb
                                - self.params.k_mAb_plasma_brain*plasma_monomer_mAb)

        d_plasma_oligomer = (self.general.k_oligomer_brain_plasma*(Visf/Vplasma)*brain_oligomer
                             - self.general.k_oligomer_plasma_brain*plasma_oligomer
                             + self.general.k_oligomer_csf_plasma*(Vcsf/Vplasma)*csf_oligomer
                             - self.params.onPP*plasma_ab*plasma_oligomer
                             + self.params.offma1*plasma_oligomer_mAb
                             - self.general.k_clear_oligomer_plasma*plasma_oligomer
                             + self.general.k_olig_inc_ext*plasma_monomer
                             - self.general.k_olig_sep_ext*plasma_oligomer)

        d_plasma_oligomer_mAb = (self.params.onPP*plasma_ab*plasma_oligomer
                                 - self.params.offma1*plasma_oligomer_mAb
                                 - self.params.plasma_clearance*plasma_oligomer_mAb
                                 + self.params.k_mAb_csf_plasma*(Vcsf/Vplasma)*csf_oligomer_mAb
                                 + self.params.k_mAb_brain_plasma*(Visf/Vplasma)*brain_oligomer_mAb
                                 - self.params.k_mAb_plasma_csf*plasma_oligomer_mAb
                                 - self.params.k_mAb_plasma_brain*plasma_oligomer_mAb)

        d_csf_monomer = (self.general.k_monomer_brain_csf*(Visf/Vcsf)*brain_monomer
                         + self.general.k_monomer_plasma_csf*(Vplasma/Vcsf)*plasma_monomer
                         - self.general.k_monomer_csf_brain*csf_monomer
                         - self.general.k_monomer_csf_plasma*csf_monomer
                         - self.params.onPP*csf_ab*csf_monomer
                         + self.params.offma0*csf_monomer_mAb
                         - self.general.k_olig_inc_ext*csf_monomer
                         + self.general.k_olig_sep_ext*csf_oligomer)

        d_csf_monomer_mAb = (self.params.onPP*csf_ab*csf_monomer
                             - self.params.offma0*csf_monomer_mAb
                             - self.params.k_mAb_csf_plasma*csf_monomer_mAb
                             + self.params.k_mAb_brain_csf*(Visf/Vcsf)*brain_monomer_mAb
                             + self.params.k_mAb_plasma_csf*(Vplasma/Vcsf)*plasma_monomer_mAb
                             - self.params.k_mAb_csf_brain*csf_monomer_mAb)

        d_csf_oligomer = (self.general.k_oligomer_brain_csf*(Visf/Vcsf)*brain_oligomer
                          - self.general.k_oligomer_csf_brain*csf_oligomer
                          - self.general.k_oligomer_csf_plasma*csf_oligomer
                          + self.general.k_oligomer_plasma_csf*(Vplasma/Vcsf)*plasma_oligomer
                          - self.params.onPP*csf_ab*csf_oligomer
                          + self.params.offma1*csf_oligomer_mAb
                          + self.general.k_olig_inc_ext*csf_monomer
                          - self.general.k_olig_sep_ext*csf_oligomer)

        d_csf_oligomer_mAb = (self.params.onPP*csf_ab*csf_oligomer
                              - self.params.offma1*csf_oligomer_mAb
                              - self.params.k_mAb_csf_plasma*csf_oligomer_mAb
                              + self.params.k_mAb_brain_csf*(Visf/Vcsf)*brain_oligomer_mAb
                              + self.params.k_mAb_plasma_csf*(Vplasma/Vcsf)*plasma_oligomer_mAb
                              - self.params.k_mAb_csf_brain*csf_oligomer_mAb)

        dYdt = [d_brain_mAb, d_brain_monomer, d_brain_oligomer, d_brain_plaque,
                d_plasma_mAb, d_plasma_monomer, d_plasma_oligomer,
                d_csf_mAb, d_csf_monomer, d_csf_oligomer,
                d_brain_monomer_mAb,  d_brain_oligomer_mAb, d_brain_plaque_mAb,
                d_plasma_monomer_mAb, d_plasma_oligomer_mAb,
                d_csf_monomer_mAb, d_csf_oligomer_mAb]

        return dYdt

    def equations_healthy(self, t, y):
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
        Vplasma = 3.5
        Vcsf = 0.15
        Visf = 0.3

        brain_ab = y[0]
        brain_monomer = y[1]
        brain_oligomer = y[2]
        brain_plaque = y[3]

        plasma_ab = y[4]
        plasma_monomer = y[5]
        plasma_oligomer = y[6]

        csf_ab = y[7]
        csf_monomer = y[8]
        csf_oligomer = y[9]

        brain_monomer_mAb = y[10]
        brain_oligomer_mAb = y[11]
        brain_plaque_mAb = y[12]

        plasma_monomer_mAb = y[13]
        plasma_oligomer_mAb = y[14]
        csf_monomer_mAb = y[15]
        csf_oligomer_mAb = y[16]

        d_plasma_mAb = (self.dosefn(self.gpm, self.dose_list, t)
                        - self.params.plasma_clearance*plasma_ab
                        + self.params.k_mAb_csf_plasma*(Vcsf/Vplasma)*csf_ab
                        - self.params.k_mAb_plasma_csf*plasma_ab
                        + self.params.k_mAb_brain_plasma*(Visf/Vplasma)*brain_ab
                        - self.params.k_mAb_plasma_brain*plasma_ab
                        - self.params.onPP*plasma_ab*plasma_monomer
                        + self.params.offma0*plasma_monomer_mAb
                        - self.params.onPP*plasma_ab*plasma_oligomer
                        + self.params.offma1*plasma_oligomer_mAb)

        d_brain_mAb = (- self.params.brain_clearance*brain_ab
                       - self.params.k_mAb_brain_plasma*brain_ab
                       - self.params.k_mAb_brain_csf*brain_ab
                       + self.params.k_mAb_csf_brain*(Vcsf/Visf)*csf_ab
                       + self.params.k_mAb_plasma_brain*(Vplasma/Visf)*plasma_ab
                       - self.params.onPP*brain_ab*brain_monomer
                       + self.params.offma0*brain_monomer_mAb
                       - self.params.onPP*brain_ab*brain_oligomer
                       + self.params.offma1*brain_oligomer_mAb
                       - self.params.onPP*brain_ab*brain_plaque
                       + self.params.offma2*brain_plaque_mAb)

        d_csf_mAb = (- self.params.k_mAb_csf_plasma*csf_ab
                     + self.params.k_mAb_brain_csf*(Visf/Vcsf)*brain_ab
                     + self.params.k_mAb_plasma_csf*(Vplasma/Vcsf)*plasma_ab
                     - self.params.k_mAb_csf_brain*csf_ab
                     - self.params.onPP*csf_ab*csf_monomer
                     + self.params.offma0*csf_monomer_mAb
                     - self.params.onPP*csf_ab*csf_oligomer
                     + self.params.offma1*csf_oligomer_mAb)

        d_brain_monomer = (self.general.k_in
                           - self.general.k_olig_inc*brain_monomer
                           + self.general.k_olig_sep*brain_oligomer
                           - self.general.k_clear_Abeta_brain*brain_monomer
                           - self.general.k_monomer_brain_plasma*brain_monomer
                           + self.general.k_monomer_plasma_brain*(Vplasma/Visf)*plasma_monomer
                           - self.general.k_monomer_brain_csf*brain_monomer
                           + self.general.k_monomer_csf_brain*(Vcsf/Visf)*csf_monomer
                           - self.params.onPP*brain_ab*brain_monomer
                           + self.params.offma0*brain_monomer_mAb)

        d_brain_monomer_mAb = (self.params.onPP*brain_ab*brain_monomer
                               - self.params.offma0*brain_monomer_mAb
                               - self.params.brain_clearance*brain_monomer_mAb
                               - self.params.k_mAb_brain_plasma*brain_monomer_mAb
                               - self.params.k_mAb_brain_csf*brain_monomer_mAb
                               + self.params.k_mAb_csf_brain*(Vcsf/Visf)*csf_monomer_mAb
                               + self.params.k_mAb_plasma_brain*(Vplasma/Visf)*plasma_monomer_mAb)

        d_brain_oligomer = (self.general.k_olig_inc*brain_monomer
                            - self.general.k_olig_sep*brain_oligomer
                            - self.general.k_plaque_inc*brain_oligomer
                            + self.general.k_plaque_sep*brain_plaque
                            - self.general.k_clear_oligomer_brain*brain_oligomer
                            - self.general.k_oligomer_brain_plasma*brain_oligomer
                            + self.general.k_oligomer_plasma_brain*(Vplasma/Visf)*plasma_oligomer
                            - self.general.k_oligomer_brain_csf*brain_oligomer
                            + self.general.k_oligomer_csf_brain*(Vcsf/Visf)*csf_oligomer
                            - self.params.onPP*brain_ab*brain_oligomer
                            + self.params.offma1*brain_oligomer_mAb)

        d_brain_oligomer_mAb = (self.params.onPP*brain_ab*brain_oligomer
                                - self.params.offma1*brain_oligomer_mAb
                                - self.params.brain_clearance*brain_oligomer_mAb
                                - self.params.k_ADCP*brain_oligomer_mAb
                                - self.params.k_mAb_brain_plasma*brain_oligomer_mAb
                                - self.params.k_mAb_brain_csf*brain_oligomer_mAb
                                + self.params.k_mAb_csf_brain*(Vcsf/Visf)*csf_oligomer_mAb
                                + self.params.k_mAb_plasma_brain*(Vplasma/Visf)*plasma_oligomer_mAb)

        d_brain_plaque = (self.general.k_plaque_inc*brain_oligomer
                          - self.general.k_plaque_sep*brain_plaque
                          - self.params.onPP*brain_ab*brain_plaque
                          + self.params.offma2*brain_plaque_mAb
                          - self.general.k_clear_plaque*brain_plaque)

        d_brain_plaque_mAb = (self.params.onPP*brain_ab*brain_plaque
                              - self.params.offma2*brain_plaque_mAb
                              - (self.params.k_ADCP)*brain_plaque_mAb)

        d_plasma_monomer = (self.general.k_peripheral_production
                            + self.general.k_monomer_brain_plasma*(Visf/Vplasma)*brain_monomer
                            - self.general.k_monomer_plasma_brain*plasma_monomer
                            + self.general.k_monomer_csf_plasma*(Vcsf/Vplasma)*csf_monomer
                            - self.params.onPP*plasma_ab*plasma_monomer
                            + self.params.offma0*plasma_monomer_mAb
                            - self.general.k_clear_Abeta_plasma*plasma_monomer
                            - self.general.k_olig_inc_ext*plasma_monomer
                            + self.general.k_olig_sep_ext*plasma_oligomer)

        d_plasma_monomer_mAb = (self.params.onPP*plasma_ab*plasma_monomer
                                - self.params.offma0*plasma_monomer_mAb
                                - self.params.plasma_clearance*plasma_monomer_mAb
                                + self.params.k_mAb_csf_plasma*(Vcsf/Vplasma)*csf_monomer_mAb
                                + self.params.k_mAb_brain_plasma*(Visf/Vplasma)*brain_monomer_mAb
                                - self.params.k_mAb_plasma_csf*plasma_monomer_mAb
                                - self.params.k_mAb_plasma_brain*plasma_monomer_mAb)

        d_plasma_oligomer = (self.general.k_oligomer_brain_plasma*(Visf/Vplasma)*brain_oligomer
                             - self.general.k_oligomer_plasma_brain*plasma_oligomer
                             + self.general.k_oligomer_csf_plasma*(Vcsf/Vplasma)*csf_oligomer
                             - self.params.onPP*plasma_ab*plasma_oligomer
                             + self.params.offma1*plasma_oligomer_mAb
                             - self.general.k_clear_oligomer_plasma*plasma_oligomer
                             + self.general.k_olig_inc_ext*plasma_monomer
                             - self.general.k_olig_sep_ext*plasma_oligomer)

        d_plasma_oligomer_mAb = (self.params.onPP*plasma_ab*plasma_oligomer
                                 - self.params.offma1*plasma_oligomer_mAb
                                 - self.params.plasma_clearance*plasma_oligomer_mAb
                                 + self.params.k_mAb_csf_plasma*(Vcsf/Vplasma)*csf_oligomer_mAb
                                 + self.params.k_mAb_brain_plasma*(Visf/Vplasma)*brain_oligomer_mAb
                                 - self.params.k_mAb_plasma_csf*plasma_oligomer_mAb
                                 - self.params.k_mAb_plasma_brain*plasma_oligomer_mAb)

        d_csf_monomer = (self.general.k_monomer_brain_csf*(Visf/Vcsf)*brain_monomer
                         + self.general.k_monomer_plasma_csf*(Vplasma/Vcsf)*plasma_monomer
                         - self.general.k_monomer_csf_brain*csf_monomer
                         - self.general.k_monomer_csf_plasma*csf_monomer
                         - self.params.onPP*csf_ab*csf_monomer
                         + self.params.offma0*csf_monomer_mAb
                         - self.general.k_olig_inc_ext*csf_monomer
                         + self.general.k_olig_sep_ext*csf_oligomer)

        d_csf_monomer_mAb = (self.params.onPP*csf_ab*csf_monomer
                             - self.params.offma0*csf_monomer_mAb
                             - self.params.k_mAb_csf_plasma*csf_monomer_mAb
                             + self.params.k_mAb_brain_csf*(Visf/Vcsf)*brain_monomer_mAb
                             + self.params.k_mAb_plasma_csf*(Vplasma/Vcsf)*plasma_monomer_mAb
                             - self.params.k_mAb_csf_brain*csf_monomer_mAb)

        d_csf_oligomer = (self.general.k_oligomer_brain_csf*(Visf/Vcsf)*brain_oligomer
                          - self.general.k_oligomer_csf_brain*csf_oligomer
                          - self.general.k_oligomer_csf_plasma*csf_oligomer
                          + self.general.k_oligomer_plasma_csf*(Vplasma/Vcsf)*plasma_oligomer
                          - self.params.onPP*csf_ab*csf_oligomer
                          + self.params.offma1*csf_oligomer_mAb
                          + self.general.k_olig_inc_ext*csf_monomer
                          - self.general.k_olig_sep_ext*csf_oligomer)

        d_csf_oligomer_mAb = (self.params.onPP*csf_ab*csf_oligomer
                              - self.params.offma1*csf_oligomer_mAb
                              - self.params.k_mAb_csf_plasma*csf_oligomer_mAb
                              + self.params.k_mAb_brain_csf*(Visf/Vcsf)*brain_oligomer_mAb
                              + self.params.k_mAb_plasma_csf*(Vplasma/Vcsf)*plasma_oligomer_mAb
                              - self.params.k_mAb_csf_brain*csf_oligomer_mAb)

        dYdt = [d_brain_mAb, d_brain_monomer, d_brain_oligomer, d_brain_plaque,
                d_plasma_mAb, d_plasma_monomer, d_plasma_oligomer,
                d_csf_mAb, d_csf_monomer, d_csf_oligomer,
                d_brain_monomer_mAb,  d_brain_oligomer_mAb, d_brain_plaque_mAb,
                d_plasma_monomer_mAb, d_plasma_oligomer_mAb,
                d_csf_monomer_mAb, d_csf_oligomer_mAb]

        return dYdt

    def dosefn(self, gpm, dose_list, t):
        """Function to convert dosing schedule into continuous function.

        Parameters
        ----------
        gpm : float
            Grams per mole of antibody
        dose_list : List
            List of dosing times in dosing schedule
        t : float
            Current time

        Returns
        ----------
        sol : float
            Input antibody concentration in nM at current time

        """

        f = 0.5
        delta = 0.1

        sol = 0
        for i, n in enumerate(dose_list):
            infusion = (((((self.dose[i] * self.bw)/1000/3)/gpm))*1e9)
            if t >= (n-0.5) and t <= (n+1.5):
                sol = ((infusion/2)/math.atan(1/delta))*(math.atan(math.sin(
                                    2*math.pi*(t-n-0.5)*f)/delta)) + infusion/2

        return sol
