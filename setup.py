from distutils.core import setup
setup(
    name = "hamnonlineng",
    packages = ["hamnonlineng"],
    version = "0.1.0a2",
    description = "Engineering Hamiltonians through Nonlinearities",
    author = "Stefan Krastanov",
    author_email = "stefan@krastanov.org",
    url = "https://github.com/hamnonlineng/hamnonlineng",
    download_url ="https://github.com/hamnonlineng/hamnonlineng/archive/master.zip",
    install_requires=['future'],
    keywords = ["physics", "quantum", "BBQ", "nonlinearity", "constraint programming", "linear programming"],
    classifiers = [
        "Programming Language :: Python",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        "Development Status :: 2 - Pre-Alpha",
        "Environment :: Other Environment",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Software Development :: Libraries :: Python Modules",
        ],
    long_description = """\
Engineering Hamiltonians through Nonlinearities
-----------------------------------------------

Following the formalism described in "Black-box superconducting circuit
quantization" http://arxiv.org/abs/1204.0587 we can expand `sin(a+b+...+h.c.)`
or `cos(a+b+...+h.c.)` or any other nonlinear `f(a+b+...+h.c.)` in order to
obtain non-linear terms in the Hamiltonian (a, b, and other letters represent
the annihilation operators of various oscillator modes). If we drive some of
the modes classically we can create any non-linear term in the effective
Hamiltonian. This module helps solve the constraint on the frequencies of the
modes (both the requirement that some of the monomials in the expansion are
resonant and that all the other terms in the expansion of the nonlinearity are
off-resonant).
"""
)
