import setuptools

import re
VERSIONFILE="phazap/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))


setuptools.setup(
    name="phazap",
    version=verstr,
    author="Jose María Ezquiaga",
    author_email="jose.ezquiaga@nbi.ku.dk",
    description="a package to identify strongly lensed gravitational wave signals rapidly",
    long_description="Gravitational wave phase reconstruction for low-latency identification of strongly lensed signals",
    packages=[
        "phazap",
    ],
    install_requires=[
        "bilby",
        "bilby_pipe",
        "pesummary",
        "getdist",
        "tqdm",
        "configargparse",
    ],
    classifiers=[
        "Programming Language :: Python :: 3.7",
    ],
    python_requires='>=3.7',
)
