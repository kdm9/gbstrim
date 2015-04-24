from setuptools import setup
import versioneer

versioneer.VCS = 'git'
versioneer.versionfile_source = '_version.py'
versioneer.versionfile_build = '_version.py'
versioneer.tag_prefix = ''
versioneer.parentdir_prefix = 'gbstrim-'

desc = """
gbstrim: Trim GBS reads of adaptors and reduce over-inflation of allele counts.
"""

install_requires = [
    "scikit-bio==0.2.3",
    "screed==0.8",
    "docopt>=0.6",
]

test_requires = [
    "coverage==3.7.1",
    "nose==1.3.6",
    "pep8==1.6.2",
]

setup(
    name="gbstrim",
    scripts=['gbstrim.py'],
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    install_requires=install_requires,
    tests_require=test_requires,
    description=desc,
    author="Kevin Murray",
    author_email="spam@kdmurray.id.au",
    url="https://github.com/kdmurray91/gbstrim",
    keywords=["bioinformatics", "sequence", "DNA", "GBS"],
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 2",
        "Development Status :: 3 - Alpha",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU General Public License v3 or later " +
            "(GPLv3+)",
    ],
)
