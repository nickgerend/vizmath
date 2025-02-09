from pathlib import Path
from setuptools import setup, find_packages

setup(
    name="vizmath",
    version="0.0.42",
    description="Visualization math toolkit.",
    long_description=(Path(__file__).parent / "README.md").read_text(),
    long_description_content_type="text/markdown",
    url="https://github.com/nickgerend/vizmath",
    author="Nick Gerend",
    author_email="nickgerend@gmail.com",
    license="MIT",
    classifiers=[
        "Programming Language :: Python :: 3.9",
        "License :: OSI Approved :: MIT License",
    ],
    packages=find_packages(exclude=['vizmath.examples*','vizmath.tests*']),
    install_requires=["numpy", "scipy", "matplotlib", "pandas"],
    include_package_data=True,
    package_data={'': ['data/*.csv']},
)