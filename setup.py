from setuptools import setup, find_packages

setup(
    name='HN-filtration',
    version="0.1.0",
    author='Marcf',
    license='MIT',
    install_requires=['numpy>=2.3.1','scipy>=1.10.0','cfractions>=2.4.1','networkx>=2.6.2','matplotlib>=3.10.7'],
    python_requires=">=3.8.0",
    packages=find_packages(),
)
