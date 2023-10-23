import setuptools

verstr = "0.0.1"

setuptools.setup(
    name="phazap",
    version=verstr,
    author="Jose María Ezquiaga",
    author_email="jose.ezquiaga@nbi.ku.dk",
    description="Gravitational wave phase reconstruction for low-latency identification of strongly lensed signals",
    packages=[
        "phazap",
    ],
    install_requires=[
        "bilby",
        "bilby_pipe",
        "pesummary",
        "getdist",
        "configargparse",
    ],
    classifiers=[
        "Programming Language :: Python :: 3.7",
    ],
    python_requires='>=3.7',
)
