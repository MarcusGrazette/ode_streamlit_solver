from sympy import symbols, Eq, solve, exp, cos, sin, Rational, latex, diff, simplify, Add, re, im

x = symbols('x')


def solve_homogeneous_ode(a, b, c):
    lambda_sym = symbols('lambda')
    aux_eq = a * lambda_sym**2 + b * lambda_sym + c
    roots = solve(aux_eq, lambda_sym)

    C1, C2 = symbols('C1 C2')
    yc_sym = None
    latex_yc = ""
    latex_roots = ""

    if len(roots) == 2:
        lambda1, lambda2 = roots
        if im(lambda1) != 0 or im(lambda2) != 0:
            alpha = re(lambda1)
            beta = im(lambda1)
            yc_sym = exp(alpha * x) * (C1 * cos(beta * x) + C2 * sin(beta * x))
            latex_yc = rf"y_c(x) = C_1 e^{{{latex(alpha)}x}} \cos({latex(beta)}x) + C_2 e^{{{latex(alpha)}x}} \sin({latex(beta)}x)"
            latex_roots = rf"Roots: $\lambda = {latex(alpha)} \pm {latex(beta)}i$ (complex conjugate)"
        elif abs(lambda1 - lambda2) < 1e-9:
            yc_sym = (C1 + C2 * x) * exp(lambda1 * x)
            latex_yc = rf"y_c(x) = (C_1 + C_2 x)e^{{{latex(lambda1)}x}}"
            latex_roots = rf"Root: $\lambda = {latex(lambda1)}$ (repeated real)"
        else:
            yc_sym = C1 * exp(lambda1 * x) + C2 * exp(lambda2 * x)
            latex_yc = rf"y_c(x) = C_1 e^{{{latex(lambda1)}x}} + C_2 e^{{{latex(lambda2)}x}}"
            latex_roots = rf"Roots: $\lambda_1 = {latex(lambda1)}$, $\lambda_2 = {latex(lambda2)}$ (distinct real)"
    elif len(roots) == 1:
        lambda1 = roots[0]
        yc_sym = (C1 + C2 * x) * exp(lambda1 * x)
        latex_yc = rf"y_c(x) = (C_1 + C_2 x)e^{{{latex(lambda1)}x}}"
        latex_roots = rf"Root: $\lambda = {latex(lambda1)}$ (repeated real)"
    else:
        raise ValueError("Unable to determine roots of the auxiliary equation.")

    return yc_sym, latex_yc, latex_roots, roots


def construct_rhs(f_type, f_params):
    if f_type == 'p':
        n = f_params['degree']
        coeffs = f_params['coeffs']
        if len(coeffs) != n + 1:
            raise ValueError("Incorrect number of polynomial coefficients.")
        f_sym = sum(Rational(coeffs[i]) * x**(n - i) for i in range(n + 1))
        yp_coeffs = symbols(f'a:{n+1}')
        yp_sym = sum(yp_coeffs[i] * x**i for i in range(n + 1))
    elif f_type == 'e':
        k = Rational(f_params['k'])
        A = Rational(f_params['A'])
        f_sym = A * exp(k * x)
        A_sym = symbols('A')
        yp_sym = A_sym * exp(k * x)
        yp_coeffs = [A_sym]
    elif f_type == 't':
        k = Rational(f_params['k'])
        a = Rational(f_params['a'])
        b = Rational(f_params['b'])
        f_sym = a * cos(k * x) + b * sin(k * x)
        A_trig, B_trig = symbols('A_trig B_trig')
        yp_sym = A_trig * cos(k * x) + B_trig * sin(k * x)
        yp_coeffs = [A_trig, B_trig]
    else:
        raise ValueError("Invalid function type.")

    return f_sym, yp_sym, yp_coeffs


def adjust_for_overlap(yp_sym, yc_sym):
    while any(simplify(term / t).is_constant() for term in Add.make_args(yp_sym) for t in Add.make_args(yc_sym)):
        yp_sym *= x
    return yp_sym


def find_particular_solution(a, b, c, yc_sym, f_type, f_params):
    f_sym, yp_sym, yp_coeffs = construct_rhs(f_type, f_params)
    yp_sym = adjust_for_overlap(yp_sym, yc_sym)

    yp_prime = diff(yp_sym, x)
    yp_double_prime = diff(yp_prime, x)

    ode_lhs = a * yp_double_prime + b * yp_prime + c * yp_sym
    equation = Eq(ode_lhs, f_sym)

    sol = solve(equation, yp_coeffs, dict=True)
    if not sol:
        raise ValueError("No solution found with assumed form.")

    sol = sol[0]
    yp_particular = yp_sym.subs(sol)
    latex_yp = f"y_{{p}}(x) = {latex(yp_particular)}"

    return yp_particular, latex_yp
