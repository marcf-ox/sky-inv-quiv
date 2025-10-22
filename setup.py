from setuptools import setup, find_packages

setup(
    name='sky-inv-quiv',
    version="0.1.1",
    author='Marcf',
    license='MIT',
    install_requires=['numpy>=2.1.1','scipy>=1.9.0','cfractions>=2.3.1','networkx>=2.5.1','matplotlib>=3.9.1'],
    python_requires=">=3.9.0",
    packages=find_packages(),
)
