from setuptools import find_packages, setup
setup(
    name=’mypythonlib’,
    packages=find_packages(include=[‘Peptide_Chef’]),
    version=’0.1.0',
    description=’My first Python library’,
    author=’TTCooper-PhD’,
    license=’MIT’,
    install_requires=[],
    setup_requires=[‘pyteomics,gzip’],
)