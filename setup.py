import os

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

# get text of README.md
current_path = os.path.dirname(os.path.realpath(__file__))
with open(os.path.join(current_path, "README.md")) as f:
    readme_text = f.read()

# get version
with open(os.path.join(current_path, "pdb4all", "__init__.py")) as f:
    for line in f:
        if line.startswith("__version__"):
            version = line.split("=")[1].strip().strip("'")
            break

setup(
    name="pdb4all",
    version=version,
    description="Convert between common protein pdb formats and names ",
    long_description=readme_text,
    long_description_content_type="text/markdown",
    url="https://github.com/boneta/pdb4all",
    author="Sergio Boneta",
    author_email="boneta@unizar.es",
    license="GPLv3",
    python_requires='>=3.6',
    classifiers=[
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
        ],
    packages=["pdb4all"],
    include_package_data=True,
    entry_points={
        "console_scripts": ["pdb4all=pdb4all.__main__:main"]
        },
    )
