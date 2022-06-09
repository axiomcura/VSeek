from setuptools import setup, find_packages

setup(
    name = 'name',
		description="simple description",
    packages = find_packages(),
		install_requires=["package1", ["package2==1.0"],
)