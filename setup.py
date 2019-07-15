# Writing the setup script: https://docs.python.org/3/distutils/setupscript.html

try:
    from setuptools import setuptools

except ImportError:
    from distutils.core import setup

config = {
    'description': 'My Project',
    'author': 'John Cooper',
    'url': 'URL to get it at.',
    'download_url': 'Where to download it.',
    'author_email': 'johnphilipcooper@gmail.com.',
    'version': '1.0',
    'install_requires': ['nose'],
    'packages': [''],
    'scripts': [],
    'project': ''
}

setup(**config)
