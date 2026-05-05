from setuptools import setup, find_packages

setup(
    name="tandem",
    version="0.2.0",
    packages=find_packages(),
    install_requires=[
        "biopython>=1.79",
        "numpy>=1.20",
    ],
    entry_points={
        "console_scripts": [
            "tandem=tandem.cli:main",
        ],
    },
    python_requires=">=3.8",
)
