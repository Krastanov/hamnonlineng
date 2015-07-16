try:
    import future
    from future.standard_library import install_aliases
    install_aliases()
except ImportError:
    pass


import fractions
import functools
import itertools
import math
import numbers
import operator


class Monomial(object):
    """A class representing a monomial of the type ``number*operators + h.c.``.

    *The h.c. part is implicit*.

    Operators have to be single letters.
    Capital letters stand for the hermitian conjugates of lower letters.

    `Monomial.n` is the numerical factor.
    `Monomial.s` is the string of operators.

    >>> Monomial(1,'a')*Monomial(1,'b')
    1.ab
    1.aB
    """
    def __new__(cls, number, symbol):
        if number == 0:
            return 0
        if symbol == '':
            return number
        symbol = Monomial._proper_sorting_string(symbol)
        normal_symbols = cls._normal_order_string(symbol)
        if len(normal_symbols) == 1:
            symbol = normal_symbols[0]
            conj_symbol = Monomial._hermitian_conjugate_string(symbol)
            if hash(conj_symbol) < hash(symbol):
                symbol = conj_symbol
            self = super(Monomial, cls).__new__(cls)
            self.n = number
            self.s = symbol
            return self
        return Polynomial([Monomial(number, _) for _ in normal_symbols])
    def __eq__(self, other):
        return isinstance(other, Monomial) and self.n == other.n and self.s == other.s
    def __hash__(self):
        return hash((self.n, self.s))
    def __add__(self, other):
        if other == self:
            return 2*self
        else:
            return Polynomial([self, other])
    def __radd__(self, other):
        return self+other
    def __mul__(self, other):
        if isinstance(other, Monomial):
            return Polynomial([Monomial(self.n*other.n, self.s+other.s),
                               Monomial(self.n*other.n, self.s+Monomial._hermitian_conjugate_string(other.s))])
        elif isinstance(other, numbers.Number):
            return Monomial(self.n*other, self.s)
        else:
            raise NotImplementedError()
    def __rmul__(self, other):
        if isinstance(other, numbers.Number):
            return Monomial(self.n*other, self.s)
        else:
            raise NotImplementedError()
    def __pow__(self, other):
        if isinstance(other, numbers.Number):
            return functools.reduce(operator.mul, [self]*other)
        else:
            raise NotImplementedError()
    def __truediv__(self, other):
        if isinstance(other, numbers.Number):
            return Monomial(self.n/other, self.s)
        else:
            raise NotImplementedError()
    def __div__(self, other):
        return self.__truediv__(other)
    def __str__(self):
        return '%s.%s'%(self.n, self.s)
    def __repr__(self):
        return str(self)
    @staticmethod
    def _normal_order_string(string):
        if not string:
            return [string]
        last_letter = string[0]
        for i, letter in enumerate(string[1:]):
            if last_letter.isupper() or letter.islower() or last_letter.lower() != letter.lower():
                last_letter = letter
                continue
            elif letter == last_letter.capitalize():
                i += 1
                return Monomial._normal_order_string(string[:i-1]+string[i]+string[i-1]+string[i+1:])\
                     + Monomial._normal_order_string(string[:i-1]+string[i+1:])
            else:
                raise RuntimeError('What!? This should not happen!')
        return [string]
    @staticmethod
    def _proper_sorting_string(string):
        return ''.join(sorted(string, key=str.lower))
    @staticmethod
    def _hermitian_conjugate_string(string):
        return Monomial._proper_sorting_string(string[::-1].swapcase())
    def dag(self):
        return Monomial(self.n, Monomial.hermitian_conjugate_string(self.s))


class Polynomial(object):
    """A class representing a sum of monomials (each monomial implicitly having its h.c.).

    `Polynomial.m` gives the list of monomials."""
    def __new__(cls, monomials):
        not_flat, flat = split_by_predicate(lambda _: isinstance(_, Polynomial), monomials)
        monomials = itertools.chain(flat, *(_.m for _ in not_flat))
        nums, not_nums = split_by_predicate(lambda _: isinstance(_, numbers.Number), monomials)
        constant = sum(nums)
        not_nums_monomials = sorted(sorted(not_nums, key=lambda _: _.s), key=lambda _: _.s.lower())
        reduced_monomials = []
        for key, group in itertools.groupby(not_nums_monomials, lambda _: _.s):
            num_factor = sum(_.n for _ in group)
            if num_factor:
                reduced_monomials.append(Monomial(num_factor, key))
        if constant:
            reduced_monomials.append(constant)
        if len(reduced_monomials) == 1:
            return reduced_monomials[0]
        elif len(reduced_monomials) == 0:
            return 0
        self = super(Polynomial, cls).__new__(cls)
        self.m = reduced_monomials
        return self
    def __add__(self, other):
        if isinstance(other, Polynomial):
            return Polynomial(self.m+other.m)
        else:
            return Polynomial(self.m+[other])
    def __radd__(self, other):
        return self+other
    def __mul__(self, other):
        if not isinstance(other, Polynomial):
            return Polynomial(_*other for _ in self.m)
        else:
            def _monomial_mul_(l, r):
                if isinstance(l, Monomial) and isinstance(r, Monomial):
                    return [Monomial(l.n*r.n, l.s+r.s),
                            Monomial(l.n*r.n, l.s+Monomial._hermitian_conjugate_string(r.s))]
                else:
                    return [l*r]
            return Polynomial(itertools.chain.from_iterable(_monomial_mul_(*_)
                                                            for _ in itertools.product(self.m, other.m)))
    def __pow__(self, other):
        if isinstance(other, numbers.Number):
            return functools.reduce(operator.mul, [self]*other)
        else:
            raise NotImplementedError()
    def __truediv__(self, other):
        if isinstance(other, numbers.Number):
            return Polynomial(_/other for _ in self.m)
        else:
            raise NotImplementedError()
    def __div__(self, other):
        return self.__truediv__(other)
    def __str__(self):
        return '\n'.join(map(str, self.m))
    def __repr__(self):
        return str(self)


def split_by_predicate(p, l):
    '''Separate an iterable in two iterators depending on predicate.'''
    t1, t2 = itertools.tee(l)
    return filter(p, t1), itertools.filterfalse(p, t2)


def operator_sum(string):
    '''Given a string of letters return the sum of operators represented by those letters.'''
    return sum((Monomial(1,_) for _ in string), 0)


def sin_terms(x, order):
    '''Given an argument return the `order`-th term in the sine's expansion.'''
    if order%2:
        if (order//2)%2:
            f = -1
        else:
            f = 1
        return x**order/(fractions.Fraction(math.factorial(order))*f)
    return 0


def drop_matching(monomials, to_drop):
    """From the list `monomials` filter out Monomials present in the list
    `to_drop`. Compare only the string part, neglect the numerical factor."""
    return itertools.filterfalse(lambda _:any(d.s==_.s for d in to_drop),
                                 monomials)


def retain_matching(monomials, to_retain):
    """The inverse of drop_matching."""
    return (_ for _ in monomials if any(d.s==_.s for d in to_retain))


def drop_definitely_offresonant(monomials):
    """From the list `monomials` filter out Monomials that contain only
    annihilation or only creation operators."""
    return (_ for _ in monomials if not (_.s.isupper() or _.s.islower()))


def drop_single_mode(monomials):
    """From the list `monomials` filter out Monomials that contain operators
    for only a single mode."""
    return itertools.filterfalse(lambda _: all(_.s[0].lower()==__ for __ in _.s[1:].lower()),
                                 monomials)


def monomial_to_constraint(monomial, letters): # XXX Assumes only Monomial instances are given to it.
    '''Given a monomial return the multiplicity of each letter, counting conjugates as negatives.'''
    constr = []
    for l in letters:
        constr.append(monomial.s.count(l)-monomial.s.count(l.upper()))
    return constr


def monomials_to_matrix(monomials, letters):
    """Run `monomial_to_constraint` on each element of the list `monomials` and
    generate the matrix of constraints."""
    import numpy as np
    return np.array([monomial_to_constraint(_, letters) for _ in monomials],dtype=int)


def solve_constraints_gecode(resonant_monomials, off_resonant_monomials, letters, maxfreq=10, elastic=0):
    '''Uses gecode to solve the constraints.'''
    import gecode
    s = gecode.space()
    variables = s.intvars(len(letters),1,maxfreq)
    s.distinct(variables)
    for m in resonant_monomials:
        if elastic:
            s.linear([-_ for _ in monomial_to_constraint(m, letters)], variables, gecode.IRT_LQ, elastic)
            s.linear(monomial_to_constraint(m, letters), variables, gecode.IRT_GQ, elastic)
        else:
            s.linear(monomial_to_constraint(m, letters), variables, gecode.IRT_EQ,0)
    for m in off_resonant_monomials:
        s.linear(monomial_to_constraint(m, letters), variables, gecode.IRT_NQ,0)
    s.branch(variables, gecode.INT_VAR_SIZE_MAX, gecode.INT_VAL_MIN)
    return (_.val(variables) for _ in s.search())


def solve_constraints_ortools(resonant_monomials, off_resonant_monomials, letters, maxfreq=10):
    '''Uses ortools to solve the constraints.'''
    from ortools.constraint_solver import pywrapcp
    from uuid import uuid4
    solver = pywrapcp.Solver(str(uuid4()))
    variables = [solver.IntVar(1, maxfreq, _) for _ in letters]
    for m in resonant_monomials:
        constr = monomial_to_constraint(m, letters)
        solver.Add(sum(c*v for c, v in zip(constr, variables)) == 0)
    for m in off_resonant_monomials:
        constr = monomial_to_constraint(m, letters)
        solver.Add(sum(c*v for c, v in zip(constr, variables)) != 0)
    for l, r in itertools.combinations(variables, 2):
        solver.Add(l!=r)
    db = solver.Phase(variables,
                      solver.INT_VAR_SIMPLE,
                      solver.INT_VALUE_SIMPLE)
    solver.NewSearch(db)
    while solver.NextSolution():
        yield [_.Value() for _ in variables]
    solver.EndSearch()


def head_and_count(iterator):
    """Print the first 10 and store the first 100000 solutions."""
    count = 0
    res = []
    for _ in iterator:
        res.append(_)
        count += 1
        if count > 100000:
            print('More than 100k solutions.')
            break
        if count < 10:
            print(_)
    else:
        print('%d solutions.'%count)
    return res


def solve_linearprog_ortools_glop(resonant_monomials, off_resonant_monomials, letters, maxfreq=10,
                                  detune=0.5):
    '''Uses ortools to solve the linear programming problem.'''
    from ortools.linear_solver import pywraplp
    from uuid import uuid4
    solver = pywraplp.Solver(str(uuid4()),
                             pywraplp.Solver.GLOP_LINEAR_PROGRAMMING,
                             #pywraplp.Solver.CLP_LINEAR_PROGRAMMING,
                            )
    variables = [solver.NumVar(0, maxfreq, _) for _ in letters]
    for m in resonant_monomials:
        constr = monomial_to_constraint(m, letters)
        solver.Add(sum(c*v for c, v in zip(constr, variables)) == 0)
    # TODO
    for m in off_resonant_monomials:
        raise NotImplementedError
    for l, r in itertools.combinations(variables, 2):
        raise NotImplementedError
        # XXX the t1-t2 trick is not working
        #t1 = solver.NumVar(0, solver.infinity(), str(uuid4()))
        #t2 = solver.NumVar(0, solver.infinity(), str(uuid4()))
        #solver.Add(l-r==t1-t2)
        #solver.Add(t1+t2>=detune)
        # XXX the binary trick is not working
        #b = solver.IntVar(0,1,str(uuid4()))
        #M = detune + maxfreq
        #solver.Add(l-r+M*b >= detune)
        #solver.Add(-l+r+M*(1-b) >= detune)
    # END TODO
    objective = solver.Objective()
    for v in variables:
        objective.SetCoefficient(v, 1)
    objective.SetMinimization()
    status = solver.Solve()
    if status:
        print('Failure with status code %d.'%status)
    print('Found a solution.')
    return [_.solution_value() for _ in variables]


def solve_linearprog_scipy(resonant_monomials, off_resonant_monomials, letters, maxfreq=10,
                           detune=0.5):
    '''Uses scipy to solve the linear programming problem.'''
    raise NotImplementedError
    #needs sparse matrices...


def solve_linearprog_pulp(resonant_monomials, off_resonant_monomials, letters, maxfreq=10,
                          detune=0.5, distinct=True):
    '''Uses ortools to solve the linear programming problem.'''
    import pulp
    from uuid import uuid4
    problem_name = str(uuid4())
    prob = pulp.LpProblem(problem_name, pulp.LpMinimize)
    variables = [pulp.LpVariable(_, 0, maxfreq, pulp.LpContinuous) for _ in letters]
    prob += sum(variables)
    print("Adding resonant constraints...")
    for m in resonant_monomials:
        constr = monomial_to_constraint(m, letters)
        expr = sum(c*v for c, v in zip(constr, variables))
        if elastic:
            prob += (expr <= elastic)
            prob += (expr >= -elastic)
        else:
            prob += (expr == 0)
    print("Adding off-resonant constraints...")
    for m in off_resonant_monomials:
        constr = monomial_to_constraint(m, letters)
        x = sum(c*v for c, v in zip(constr, variables))
        b = pulp.LpVariable(str(uuid4()), 0, 1, pulp.LpBinary)
        M = detune + maxfreq*len(m.s)
        prob += (x+M*b >= detune)
        prob += (-x-M*b >= detune-M)
    if distinct:
        print("Adding distinct constraints...")
        for l, r in itertools.combinations(variables, 2):
            b = pulp.LpVariable(str(uuid4()), 0, 1, pulp.LpBinary)
            M = detune + maxfreq
            prob += (l-r+M*b >= detune)
            prob += (-l+r-M*b >= detune-M)
    else:
        print("Skipping distinct contraints.")
    print("Dumping LP file...")
    prob.writeLP("%s.lp"%problem_name)
    print("Solving...")
    prob.solve()
    print("Status: %d - %s"%(prob.status,pulp.LpStatus[prob.status]))
    return [(_.name, _.value()) for _ in variables]


def detuning_of_monomial(solution, monomials, letters):
    """For the frequencies given in `solution` calculate the detuning of the
    Monomials given in the list `monomials`."""
    import numpy as np
    mat = monomials_to_matrix(monomials, letters)
    det = np.abs(mat.dot(solution))
    return det


def filter_resonant(solution, monomials, letters):
    """Return only the resonant Monomials from the list `monomials` given the
    frequencies in `solution`."""
    import numpy as np
    det = detuning_of_monomial(solution, monomials, letters)
    return [monomials[i] for i in np.argwhere(det==0)]
