

from setuptools import setup, find_packages


name = 'gethog3'
__version__ = None
with open('{:s}/__init__.py'.format(name), 'rt') as fp:
    for line in fp:
        if line.startswith('__version__'):
            exec(line.rstrip())

# TODO
requirements = []

desc = 'GetHOG3 - blabla'

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
            "infer-roothogs=gethog3.infer_roothog:infer_roothogs",
            #"infer-folder=gethogs.infer_folder:infer_folder",
        ]
    },
)
