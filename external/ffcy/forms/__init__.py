from ffcy.reference import helmholtz_action, laplace_action, mass_action
from forms.helmholtz import helmholtz_form
from forms.laplace import laplace_form
from forms.mass import mass_form

# Each form with the NumPy reference for its action.
FORMS = {
    "mass": (mass_form, mass_action),
    "laplace": (laplace_form, laplace_action),
    "helmholtz": (helmholtz_form, helmholtz_action),
}
