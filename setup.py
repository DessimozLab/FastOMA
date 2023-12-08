

from setuptools import setup, find_packages

name = 'FastOMA'
__version__ = None
with open('{:s}/__init__.py'.format(name), 'rt') as fp:
    for line in fp:
        if line.startswith('__version__'):
            exec(line.rstrip())

# TODO
requirements = ['biopython', 'ete3', 'omamer>=2.0.0', 'nextflow', 'pyparsing' , 'DendroPy', 'future', 'lxml','pyham', 'networkx']

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
            "fastoma-check-input=FastOMA.check_input:fastoma_check_input",
            "fastoma-infer-roothogs=FastOMA.infer_roothogs:fastoma_infer_roothogs",
            "fastoma-batch-roothogs=FastOMA.batch_roothogs:fastoma_batch_roothogs",
            "fastoma-infer-subhogs=FastOMA.infer_subhogs:fastoma_infer_subhogs",
            "fastoma-collect-subhogs=FastOMA.collect_subhogs:fastoma_collect_subhogs",
        ]
    },
)
