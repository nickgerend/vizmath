import pathlib
import setuptools as st
from setuptools import setup

HERE = pathlib.Path(__file__).parent
README = (HERE / "README.md").read_text()
setup(
    name="vizmath",
    version="0.0.14",
    description="Visualization math toolkit.",
    long_description="Welcome to vizmath! Please visit [GitHub](https://github.com/nickgerend/vizmath/blob/main/README.md)",
    long_description_content_type="text/markdown",
    url="https://github.com/nickgerend/vizmath",
    author="Nick Gerend",
    author_email="nickgerend@gmail.com",
    license="MIT",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.9",
    ],
    packages=st.find_namespace_packages(),
    install_requires=["numpy", "scipy", "matplotlib", "pandas"],
    include_package_data=True,
    package_data={'': ['data/*.csv']},
)