from setuptools import setup, find_packages

setup(
    name="cmd_fusion",
    version="0.1.0",
    description="Fusion results comparison dashboard.",
    url="https://github.com/TDMedina/CMD_Fusion",
    author="T.D. Medina",
    author_email="tylerdanmedina@gmail.com",
    packages=find_packages(),
    zip_safe=False,
    python_requires=">=3.12",
    install_requires=[
        "dash==2.18.2",
        "dash-bootstrap-components==1.7.1",
        "liftover==1.3.2",
        "numpy==2.2.3",
        "pandas==2.2.3",
        "plotly==6.0.0",
        "pysam==0.23.0",
        "PyYAML==6.0.2",
        "tqdm==4.67.1"
        ],
    entry_points={
        "console_scripts": [
            "FusionDash = cmd_fusion.dashboard.__main__:main",
            ]
        },
    include_package_data=True
    )
