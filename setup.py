import pathlib
from setuptools import setup, find_packages

HERE = pathlib.Path(__file__).parent
README = (HERE / "README.md").read_text()
setup(
    name="vizmath",
    version="0.0.29",
    description="Visualization math toolkit.",
    long_description="Welcome to vizmath! Please visit [GitHub](https://github.com/nickgerend/vizmath/blob/main/README.md)",
    long_description_content_type="text/markdown",
    url="https://github.com/nickgerend/vizmath",
    author="Nick Gerend",
    author_email="nickgerend@gmail.com",
    license="Dual License: Non-Commercial Use License and Commercial Use License, see LICENSE-NC and LICENSE-COM for terms",
    classifiers=[
        "Programming Language :: Python :: 3.9",
    ],
    packages=find_packages(exclude=['vizmath.examples*','vizmath.tests*']),
    install_requires=["numpy", "scipy", "matplotlib", "pandas"],
    include_package_data=True,
    package_data={'': ['data/*.csv']},
)