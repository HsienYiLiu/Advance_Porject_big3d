# This dockerfile uses the ubuntu image
# VERSION 2 - EDITION 1
# Author: Hsien-Yi Liu
# Command format: Instruction [arguments / command] ..

# From fist command for importing image
FROM gcc:4.9
COPY . /home/liu/advance_porject_big3d
WORKDIR /home/liu/advance_porject_big3d
#RUN gcc -o myapp readFile.cpp
RUN gcc readFile.cpp -lstdc++ -o myapp
CMD ["./myapp"]
