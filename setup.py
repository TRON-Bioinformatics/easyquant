from setuptools import find_packages, setup
import easy_quant
from easy_quant.version import version

#VERSION = version


# parses requirements from file
with open("requirements.txt") as f:
    required = f.read().splitlines()

with open("README.md", "r", encoding="utf-8") as f:
    long_description = f.read()

# Build the Python package
setup(
    name='easy_quant',
    version=version,
    packages=find_packages(exclude=["legacy"]),
    entry_points={
        'console_scripts': [
            'easy_quant=easy_quant.command_line:easy_quant_cli',
        ],
    },
    author="TRON - Translational Oncology at the University Medical Center of the Johannes Gutenberg University Mainz"
    "- Computational Medicine group",
    author_email='patrick.sorn@tron-mainz.de',
    description='Quantification of reads at defined positions to verify custom input sequences. Given a gene fusion or splicing junction of interest, this tool can quantify RNA-seq reads supporting the breakpoint (or splice junction) by quantifying reads that map to the breakpoint (junction reads) and read pairs that span the breakpoint (spanning pairs).',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/tron-bioinformatics/easyquant",
    requires=[],
    install_requires=required,
    classifiers=[
        'Development Status :: 4 - Beta',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
        'Intended Audience :: Healthcare Industry',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3 :: Only',
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix"
      ],
    python_requires='>=3.7',
    license='MIT'
)
