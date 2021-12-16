from setuptools import setup, find_packages
import sys
import versioneer

try:
    import pandas_pareto.pareto
except:
    sys.exit("First install pandas_pareto:\n\ncd pandas_pareto\npython setup.py"
             " install\n")

setup (
    name='gmak',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'gmak=gmak.scripts.RunFromInput:main',
        ]
    },
    install_requires=['numpy', 'scikit-learn', 'alchemlyb', 'pymbar'],
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
)
