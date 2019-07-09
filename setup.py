import setuptools

requirements = [
    'numpy',
    'pandas',
    'scipy',
    'click',
]
setuptools.setup(
    name="MD2",
    version="0.1.0",
    url="https://github.com/dcdanko/MD2",
    author="Chandrima Bhattacharya",
    author_email="chb4004@med.cornell.edu",
    description="Codes for Microbe Directory 2.0 and above",
    packages=setuptools.find_packages(),
    package_dir={'MD2': 'MD2'},
    entry_points={
        'console_scripts': [
            'MD2=MD2.cli:main',
        ]
    },
    install_requires=requirements,
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
	package_data={'MD2': [
        'ncbi_tree/*.dmp.gz',
    ]},
)
