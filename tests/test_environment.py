"""Step 1 check: the Python environment has what the project needs, and the package imports."""


def test_scientific_stack_imports():
    import matplotlib
    import numpy
    import scipy
    import sympy

    assert numpy.__version__ and scipy.__version__ and matplotlib.__version__ and sympy.__version__


def test_project_package_imports():
    import fishtail

    assert "PORT" in fishtail.__doc__  # sign convention is documented at the top of the package
