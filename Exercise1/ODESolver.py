#Code for solving ODE in python from ode book
#Reference:
#https://sundnes.github.io/solving_odes_in_python/ode_book.pdf 
#pag 35


#Importing libraries
import numpy as np
from scipy.optimize import root

#Defining class ODESolver
class ODESolver:
    
    def __init__(self, f):
        # Wrap user’s f in a new function that always
        # converts list/tuple to array (or let array be array)
        self.f = lambda t, u: np.asarray(f(t, u), float)

    def set_initial_condition(self, u0):
        if isinstance(u0, (float, int)):  # scalar ODE
            self.neq = 1  # no of equations
            u0 = float(u0)
        else:  # system of ODEs
            u0 = np.asarray(u0)
            self.neq = u0.size  # no of equations
        self.u0 = u0

    def debug(self, activate):
        if activate == True:
            self.debug = True

    def solve(self, t_span, N):
        """Compute solution for
        t_span[0] <= t <= t_span[1],
        using N steps."""
        t0, T = t_span
        self.dt = T/N
        self.t = np.zeros(N+1)  # N steps ~ N+1 time points
        if self.neq == 1:
            self.u = np.zeros(N+1)
        else:
            self.u = np.zeros((N+1, self.neq))
        self.t[0] = t0
        self.u[0] = self.u0
        for n in range(N):
            self.n = n
            self.t[n+1] = self.t[n] + self.dt
            self.u[n+1] = self.advance()
            if self.debug == True:
                print(self.u[n+1])
        return self.t, self.u


class ExplicitMidpoint(ODESolver):
    def advance(self):
        u, f, n, t = self.u, self.f, self.n, self.t
        dt = self.dt
        dt2 = dt/2.0
        k1 = f(t[n], u[n])
        k2 = f(t[n] + dt2, u[n] + dt2*k1)
        unew = u[n] + dt*k2
        return unew


class RungeKutta4(ODESolver):
    def advance(self):
        u, f, n, t = self.u, self.f, self.n, self.t
        dt = self.dt
        dt2 = dt/2.0
        k1 = f(t[n], u[n],)
        k2 = f(t[n] + dt2, u[n] + dt2*k1, )
        k3 = f(t[n] + dt2, u[n] + dt2*k2, )
        k4 = f(t[n] + dt, u[n] + dt*k3, )
        unew = u[n] + (dt/6.0)*(k1 + 2*k2 + 2*k3 + k4)
        return unew


class ImplicitRK(ODESolver):
    def solve_stages(self):
        u, f, n, t = self.u, self.f, self.n, self.t
        neq = self.neq
        s = self.stages
        k0 = f(t[n], u[n])
        k0 = np.hstack([k0 for i in range(s)])
        sol = root(self.stage_eq, k0)
        return np.split(sol.x, s)

    def stage_eq(self, k_all):
        a, c = self.a, self.c
        s, neq = self.stages, self.neq
        u, f, n, t = self.u, self.f, self.n, self.t
        dt = self.dt
        res = np.zeros_like(k_all)
        k = np.split(k_all, s)
        for i in range(s):
            fi = f(t[n]+c[i]*dt, u[n]+dt*sum([a[i, j]*k[j] for j in range(s)]))
            res[i*neq:(i+1)*neq] = k[i] - fi
        return res

    def advance(self):
        b = self.b
        u, n, t = self.u, self.n, self.t
        dt = self.dt
        k = self.solve_stages()
        return u[n]+dt*sum(b_*k_ for b_, k_ in zip(b, k))


class ImpMidpoint(ImplicitRK):
    def __init__(self, f):
        super().__init__(f)
        self.stages = 1
        self.a = np.array([[1/2]])
        self.c = np.array([1/2])
        self.b = np.array([1])


class Radau2(ImplicitRK):
    def __init__(self, f):
        super().__init__(f)
        self.stages = 2
        self.a = np.array([[5/12, -1/12], [3/4, 1/4]])
        self.c = np.array([1/3, 1])
        self.b = np.array([3/4, 1/4])


class Radau3(ImplicitRK):
    def __init__(self, f):
        super().__init__(f)
        self.stages = 3
        sq6 = np.sqrt(6)
        self.a = np.array([[(88-7*sq6)/360, (296-169*sq6)/1800, (-2+3*sq6)/(225)],
                           [(296+169*sq6)/1800, (88+7*sq6)/360, (-2-3*sq6)/(225)], [(16-sq6)/36, (16+sq6)/36, 1/9]])
        self.c = np.array([(4-sq6)/10, (4+sq6)/10, 1])
        self.b = np.array([(16-sq6)/36, (16+sq6)/36, 1/9])