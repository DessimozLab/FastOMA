

from setuptools import setup, find_packages

name = 'FastOMA'
__version__ = None
with open('{:s}/__init__.py'.format(name), 'rt') as fp:
    for line in fp:
        if line.startswith('__version__'):
            exec(line.rstrip())

# TODO
requirements = ['biopython', 'ete3', 'omamer', 'nextflow', 'pyparsing' , 'DendroPy', 'future', 'lxml']

desc = 'FastOM - a package to infer orthology information '

setup(
    name=name,
    version=__version__,
    author='',
    email='',
    url='https://github.com/DessimozLab/gethog3',
    description=desc,
    packages=find_packages(),
    install_requires=requirements,
    python_requires=">=3.6",
    license='MIT',
    entry_points={
        'console_scripts': [
            "check-input-fastoma=FastOMA.check_input_fastoma:check_input_fastoma",
            "infer-roothogs=FastOMA.infer_roothogs:infer_roothogs",
            "batch-roothogs=FastOMA.batch_roothogs:batch_roothogs",
            "infer-subhogs=FastOMA.infer_subhogs:infer_subhogs",
            "collect-subhogs=FastOMA.collect_subhogs:collect_subhogs",
        ]
    },
)
