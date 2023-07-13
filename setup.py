import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name='EEWsimpy',
    version='0.0.2',
    author='Fred Massin',
    author_email='fmassin@ethz.ch',
    description='Testing installation of Package',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/ATTAC-EEW-in-Central-America/EEWsimpy',
    project_urls = {
        "Bug Tracker": "https://github.com/ATTAC-EEW-in-Central-America/EEWsimpy/issues"
    },
    license='AGPL-3.0',
    packages=['eewsimpy'],
    install_requires=['obspy','shapely','cartopy']
)