#!/usr/bin/env python3
from setuptools import setup

extras={
        'plotting': [
            'matplotlib'
        ],
        'bio': [
            'cameo'
        ],
        'lattice': [
            'scipy'
        ]
    }

setup(name='fea',
    version='1.00',
    description='A package for efficiently finding the n-dimensional reduced solution space from an m-dimensional problem',
    author='Matthew Long',
    license='MIT',
    packages=['fea'],
    install_requires=['numpy','pandas','optlang','networkx<2.0','sortedcontainers','glpk'],
    extras_requires=extras
)