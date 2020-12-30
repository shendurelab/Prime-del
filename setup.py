import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="primedel",
    version="1.0",
    author="Wei Chen",
    author_email="wchen108@uw.edu",
    description="A package for primedel",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/shendurelab/Prime-del/main/",
    packages=['primedel'],
    package_data={'primedel': ['indel_ratio.npz']},
    install_requires=['numpy','pandas','regex','xlrd'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    include_package_data=True,zip_safe=False)