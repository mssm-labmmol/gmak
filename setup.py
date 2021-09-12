from setuptools import setup, find_packages

setup (
    name='gmak',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'gmak=gmak.scripts.RunFromInput:main',
        ]
    },
    install_requires=['numpy', 'scikit-learn', 'alchemlyb', 'pymbar']
)
