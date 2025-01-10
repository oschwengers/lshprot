
from os import path
from setuptools import setup, find_packages
import lshprot


# Get the long description from the README file
setup_dir = path.abspath(path.dirname(__file__))
with open(path.join(setup_dir, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()


setup(
    name='lshprot',
    version=lshprot.__version__,
    description='lshProt: rapid identification of near-identical protein sequences',
    keywords=['bioinformatics', 'proteins', 'alignments', 'blast'],
    long_description=long_description,
    long_description_content_type='text/markdown',
    license='GPLv3',
    author='Oliver Schwengers',
    author_email='oliver.schwengers@computational.bio.uni-giessen.de',
    url='https://github.com/oschwengers/lshprot',
    packages=find_packages(include=['lshprot']),
    python_requires='>=3.8, <3.12',
    include_package_data=False,
    zip_safe=False,
    install_requires=[
        'biopython >= 1.78',
        'xopen >= 1.5.0',
        'datasketch >= 1.6.5',
        'pyopal >= 0.6.1',
        'alive-progress >= 3.0.1'
    ],
    entry_points={
        'console_scripts': [
            'lshprot=lshprot.main:main'
        ]
    },
    classifiers=[
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Development Status :: 5 - Production/Stable',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Natural Language :: English'
    ],
    project_urls={
        'Source': 'https://github.com/oschwengers/lshprot',
        'Bug Reports': 'https://github.com/oschwengers/lshprot/issues',
        'CI': 'https://github.com/oschwengers/lshprot/actions'
    },
)
