import setuptools

requirements = [
    'numpy',
    'pandas',
    'scipy',
    'click',
]
setuptools.setup(
    name="microbe_directory",
    version="0.1.0",
    url="https://github.com/dcdanko/MD2",
    author="Chandrima Bhattacharya",
    author_email="chb4004@med.cornell.edu",
    description="Codes for Microbe Directory 2.0 and above",
    packages=setuptools.find_packages(),
    package_dir={'microbe_directory': 'microbe_directory'},
    entry_points={
        'console_scripts': [
            'microbe_directory=microbe_directory.cli:main',
        ]
    },
    install_requires=requirements,
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    package_data={
        'microbe_directory': [
            'ncbi_tree/*.dmp.gz',
            'stored_final_tables/*.csv.gz',
        ]
    },
)
