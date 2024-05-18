from setuptools import setup, find_packages

# Lire le fichier __version__.py
version = {}
with open("OUREA-EasyCode/__version__.py") as fp:
    exec(fp.read(), version)

setup(
    name="OUREA-EasyCode",
    version=version['__version__'],
    author="Anthony PARISOT",
    author_email="ourea.engineering@gmail.com",
    description="Le but est de créer un catalogue de fonction permettant une utilisation rapide pour la réalisation de note de calcul personnalisée.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/Icelone73/OUREA-EasyCode",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.12',
)
