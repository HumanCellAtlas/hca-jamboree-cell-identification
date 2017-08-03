FROM ubuntu:16.04
# this designates who is responsible for maintaining this
LABEL maintainer="your.name@your.institution.org"
# add your code from the working directory (also called context)
COPY script.R your.example.R
# COPY HelloWorld.ipynb HelloWorld.ipynb
# install your dependencies
RUN apt-get -m update && apt-get install -y tree
# see extended example at https://github.com/denis-yuen/icy-blackberry
CMD ["uname", "-a"]
