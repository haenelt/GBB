import setuptools

# version from textfile?????????????????????????????????????????????????????
# install_requires=['peppercorn'],

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="gbb", # Replace with your own username
    version="0.1.0",
    author="Daniel Haenelt",
    author_email="daniel.haenelt@gmail.com",
    description="Gradient-based boundary (GBB) surface refinement",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/haenelt/GBB",
    licesen="LICENSE", # brauche ich das?????????????????????????????????????
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        ],
    python_requires='>=3',
    )




#      scripts=['bin/script1','bin/script2'],