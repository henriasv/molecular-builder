import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="molecular_builder", # Replace with your own username
    version="0.1.0",
    author="Henrik Andersen Sveinsson",
    author_email="henriasv@fys.uio.no",
    description="Package for building moleular systems",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/henriasv/molecular-builder",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)