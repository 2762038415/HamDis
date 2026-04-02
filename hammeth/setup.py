from setuptools import setup, find_packages

setup(
    name="hammeth",
    version="0.1.0",
    description="Hamming distance based methylation analysis toolkit",
    author="JYP",
    packages=find_packages(),
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "hammeth=hammeth.cli:main",
        ],
    },
    install_requires=[
        "pandas",
    ],
)
