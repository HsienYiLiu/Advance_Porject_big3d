# This dockerfile uses the ubuntu image
# VERSION 2 - EDITION 1
# Author: Hsien-Yi Liu
# Command format: Instruction [arguments / command] ..

# From fist command for importing image
FROM ubuntu

# Name of Maintainer
MAINTAINER Hsien-Yi Liu

# Updating command
RUN echo "deb http://archive.ubuntu.com/ubuntu/ raring main universe" >> /etc/apt/sources.list
RUN ../shane/hello
#RUN apt-get update && apt-get install -y nginx
#RUN echo "\ndaemon off;" >> /etc/nginx/nginx.conf

# Command for creating new container
#CMD /usr/sbin/nginx
