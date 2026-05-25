import scipy


class Solution:
    """Class for obtaining numeric solutions of models.

    """
    def __init__(self, model, t_0, t_end, step_size):
        """Constructor Method.

        Parameters
        ----------
        model : class
            Class containing equations() function comprising model ODEs
        t_0 : Int
            Start time for simulation
        t_end : Int
            Finish time for simulation
        step_size : Int
            Interval between start and end time at which to return
            simulation results

        """
        self.t_start = t_0
        self.t_end = t_end
        self.t_eval = list(range(t_0, t_end, step_size))
        self.model = model

        if self.model.type == 'no_ab':
            self.equations = model.equations
            self.y0 = [1.05, 0, 14.6, 0, 1300, 0,
                       6.936e-3, 0, 1.092e-3, 0,
                       0.15544, 0, 6.9138e-3, 0]
            self.max = 1

        elif self.model.type == 'healthy':
            self.equations = model.equations
            self.y0 = [1.6104, 0, 1.0701, 0, 0, 0,
                       2.5357e-3, 0, 1.72238e-3, 0,
                       0.123324, 0, 0.1147798, 0]
            self.max = 1

        elif self.model.type == 'one_ab':
            self.equations = model.equations
            self.y0 = [0, 1.05, 14.6, 1300,
                       0,  6.936e-3, 1.092e-3,
                       0, 0.15544, 6.9138e-3,
                       0, 0, 0,
                       0, 0,
                       0, 0]
            self.max = 1

        elif self.model.type == 'one_ab_healthy':
            self.equations = model.equations_healthy
            self.y0 = [0, 1.6104, 1.0701, 0,
                       0,  2.5357e-3, 1.72238e-3,
                       0, 0.123324, 0.1147798,
                       0, 0, 0,
                       0, 0,
                       0, 0]
            self.max = 1

    def solve(self):
        """Numeric solver.

        Returns
        ----------
        solution : List
            List containing solution to ODEs

        """
        solution = scipy.integrate.solve_ivp(fun=lambda t, y: self.equations(t, y),
                                             t_span=[self.t_eval[0],
                                                     self.t_eval[-1]],
                                             y0=self.y0,
                                             t_eval=self.t_eval,
                                             max_step=self.max,
                                             method='LSODA')

        return solution
