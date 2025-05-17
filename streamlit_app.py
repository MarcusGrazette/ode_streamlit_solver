import streamlit as st
from ode_solver import solve_homogeneous_ode, find_particular_solution

st.set_page_config(page_title="Second-Order ODE Solver", layout="centered")

st.title("🧮 Second-Order Linear ODE Solver")

st.markdown("### Equation format:")
st.latex(r"a y'' + b y' + c y = f(x)")

# --- Coefficients ---
st.sidebar.header("Equation Coefficients")
a = st.sidebar.number_input("Coefficient a", value=1.0)
b = st.sidebar.number_input("Coefficient b", value=0.0)
c = st.sidebar.number_input("Coefficient c", value=0.0)

# --- RHS Type ---
f_type = st.sidebar.selectbox("Choose RHS function type", ['Polynomial', 'Exponential', 'Trigonometric'])

# --- RHS Parameters ---
f_params = {}

if f_type == 'Polynomial':
    deg = st.sidebar.slider("Polynomial Degree", 0, 5, 1)
    coeffs = []
    for i in range(deg + 1):
        coeff = st.sidebar.number_input(f"Coefficient of x^{deg - i}", value=0.0)
        coeffs.append(coeff)
    f_params = {"degree": deg, "coeffs": coeffs}
    internal_type = 'p'

elif f_type == 'Exponential':
    A = st.sidebar.number_input("Amplitude (A)", value=1.0)
    k = st.sidebar.number_input("Exponent coefficient (k)", value=1.0)
    f_params = {"A": A, "k": k}
    internal_type = 'e'

elif f_type == 'Trigonometric':
    a = st.sidebar.number_input("Coefficient of cos(kx)", value=1.0)
    b = st.sidebar.number_input("Coefficient of sin(kx)", value=1.0)
    k = st.sidebar.number_input("k (frequency)", value=1.0)
    f_params = {"a": a, "b": b, "k": k}
    internal_type = 't'

# --- Solve ODE ---
if st.button("Solve ODE"):
    try:
        yc, latex_yc, latex_roots, roots = solve_homogeneous_ode(a, b, c)
        yp, latex_yp = find_particular_solution(a, b, c, yc, internal_type, f_params)
        y_general = yc + yp

        st.markdown("### ✅ Homogeneous Solution")
        st.latex(latex_roots)
        st.latex(latex_yc)

        st.markdown("### ✅ Particular Solution")
        st.latex(latex_yp)

        st.markdown("### ✅ General Solution")
        st.latex(f"y(x) = {latex(y_general)}")

    except Exception as e:
        st.error(f"An error occurred: {e}")
