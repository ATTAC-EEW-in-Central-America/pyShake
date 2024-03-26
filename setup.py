import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name='pyShake',
    version='1.0.1',
    author='Fred Massin',
    author_email='fmassin@ethz.ch',
    description='Testing installation of Package',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/ATTAC-EEW-in-Central-America/pyShake',
    project_urls = {
        "Bug Tracker": "https://github.com/ATTAC-EEW-in-Central-America/pyShake/issues"
    },
    license='AGPL-3.0',
    packages=['pyShake'],
    install_requires=['obspy','shapely','cartopy']
)