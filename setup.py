from setuptools import find_packages, setup

config = {
    'name': 'Acoustipy',
    'version': '0.1.0',
    'description': 'Tools for characterizing porous materials and their acoustic properties',
    'packages': find_packages(),
    'install_requires':[
        'matplotlib==3.5.3',
        'numpy==1.23.2',
        'pandas==1.5.3',
        'scipy==1.9.0',
    ]
}

setup(**config)