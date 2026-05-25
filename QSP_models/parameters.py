

class NoAbParameters:
    """Class defining the default parameters for use in class:NoAbModel.

    """
    def __init__(self):
        """Constructor Method.

        Establishes baseline parameters for disease model.

        """
        self.k_in = 10**(-0.82444)
        self.k_peripheral_production = 10**(-3.5429)
        self.k_olig_inc = 10**(-2.0571)
        self.k_olig_sep = 10**(-3.6848)
        self.k_olig_inc_ext = 10**(-1.3891)
        self.k_olig_sep_ext = 10**(-2.6005)
        self.k_plaque_inc = 10**(-3.9241)
        self.k_plaque_sep = 10**(-6.8431)
        self.k_clear_Abeta_plasma = 10**(-0.042913)
        self.k_clear_oligomer_plasma = 10**(-0.025126)
        self.k_clear_Abeta_brain = 10**(-1.1315)
        self.k_clear_oligomer_brain = 10**(-4.6748)
        self.k_monomer_plasma_brain = 10**(-9.3888)
        self.k_oligomer_plasma_brain = 10**(-8.4999)
        self.k_monomer_brain_plasma = 10**(-1.2383)
        self.k_oligomer_brain_plasma = 10**(-3.5286)
        self.k_monomer_brain_csf = 10**(-2.3059)
        self.k_oligomer_brain_csf = 10**(-6.0941)
        self.k_monomer_csf_brain = 10**(-9.3066)
        self.k_oligomer_csf_brain = 10**(-11.502)
        self.k_monomer_csf_plasma = 10**(-1.5821)
        self.k_oligomer_csf_plasma = 10**(-0.036815)
        self.k_clear_plaque = 10**(-7.2285)
        self.k_clear_Abeta_csf = 0
        self.k_clear_oligomer_csf = 0
        self.k_monomer_plasma_csf = 0
        self.k_oligomer_plasma_csf = 0


class Lecanemab:
    """Class defining the default parameters for use in class:OneAbModel.

    """
    def __init__(self):
        """Constructor Method.

        Establishes baseline PK and PK parameters.

        """
        # antibody specific
        self.onPP = 0.001
        self.offma0 = 2300*self.onPP
        self.offma1 = 67.3*self.onPP
        self.offma2 = 1.8*self.onPP

        # general parameters
        self.k_mAb_plasma_brain = 10**(-6)
        self.k_mAb_brain_csf = 10**(-0.1)
        self.k_mAb_csf_brain = 10**(-0.49683)
        self.k_mAb_csf_plasma = 10**(-2.4334)
        self.k_mAb_plasma_csf = 10**(-8.0)
        self.plasma_clearance = 10**(-2.2989)
        self.k_mAb_brain_plasma = 10**(-5.7113)
        self.brain_clearance = 0
        self.k_ADCP = 10**(-1.1958)
        self.scale = 10.52631579
