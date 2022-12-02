from setuptools import find_packages, setup


VERSION='0.1.0'


setup(
    name="Peptide_Chef",
    version=VERSION,
    author="Tyler T. Cooper",
    author_email="tcoope2@gmail.com"
    description="An all-in-one Python Libary for Proteomic Data Analyses"
    long_description=open('README.txt').read()+ '\n\n' + open('CHANGELOG.txt').read(),
    long_description_content_type="text/markdown",
    url="https://github.com/TTCooper-PhD",
    project_urls={
        "Bug Tracker": "https://github.com/TTCooper-PhD",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        'Intended Audience :: Education'
    ],
    keywords=['Proteomics',"Python"."DataAnalysis"],
    python_requires=">=3.6",
    install_requires=["json","numpy","pandas","seaborn","matplotlib","itertools","re","gzip","urllib"]
    packages=find_packages(include=["Peptide_Chef"]),
    description="Peptide Chef: A Python-based Tool for Proteomic Analyses and Datavisulization",
    license="MIT",
)