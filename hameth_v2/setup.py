from pathlib import Path

from setuptools import find_packages, setup


ROOT = Path(__file__).resolve().parent

setup(
    name="hameth-v2",
    version="2.0.0",
    description="Depth-balanced read-pattern Hamdist analysis",
    long_description=(ROOT / "README.md").read_text(encoding="utf-8"),
    long_description_content_type="text/markdown",
    url="https://github.com/2762038415/HamDis/tree/main/hameth_v2",
    maintainer="2762038415",
    maintainer_email="2762038415@users.noreply.github.com",
    packages=find_packages(),
    include_package_data=True,
    python_requires=">=3.7",
    install_requires=["numpy", "pandas", "pysam", "tqdm"],
    entry_points={
        "console_scripts": [
            "hameth_v2=hameth_v2.cli:main",
            "hameth-v2=hameth_v2.cli:main",
        ]
    },
)
