"""Install configuration for brainshuttle_QSP."""

from pathlib import Path

from setuptools import find_packages, setup


def read_requirements(filename: str = "requirements.txt") -> list[str]:
    """Read non-comment lines from requirements.txt for install_requires."""
    path = Path(__file__).parent / filename
    requirements = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        requirements.append(line)
    return requirements


setup(
    name="brainshuttle_QSP",
    version="0.1.0",
    description="Quantitative systems pharmacology models for anti-amyloid therapies",
    license="MIT",
    maintainer="Lara Herriott",
    maintainer_email="herriott@maths.ox.ac.uk",
    python_requires=">=3.10",
    packages=find_packages(include=("QSP_models", "QSP_models.*")),
    install_requires=read_requirements(),
)
