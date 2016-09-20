import os
from setuptools import setup

setup(
    name='Simple Example',
    version='0.0.1',
    author='Cannot believe I wrote this',
    author_email='dont_send_me_email@gmail.com',
    packages=['lgn_simulator'],
    license='Use as you wish. No guarantees whatsoever.',
    install_requires=[''],
    classifiers=['Development Status :: 3 - Alpha'],
    package_data={'lgn_simulator': ['*.so']},
)


