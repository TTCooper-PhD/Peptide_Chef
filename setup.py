from setuptools import find_packages, setup

setup(
    name="Peptide_Chef",
    version="0.0.7",
    author="Tyler T. Cooper",
    author_email="tcoope2@uwo.ca",
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
    keywords='Proteomics',
    python_requires=">=3.6",
    packages=find_packages(include=["Peptide_Chef"]),
    description="Peptide Chef: A Python-based Tool for Proteomic Analyses and Datavisulization",
    license="MIT",
)