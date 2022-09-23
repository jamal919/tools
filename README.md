# SIROCCO tools in Docker
The SIROCCO tools are a software suite of executables and libraries for pre and post-processings in ocean modelling (configuration set-up and data/simulations analyses). It is aimed to accept and be compliant with most Earth science formats.

It is a part of the French national COMODO (community developments for ocean modeling) initiative.

SIROCCO is funded by INSU and Observatoire Midi-Pyrénées/Université Paul Sabatier and receives project support from CNES, SHOM, IFREMER and ANR.

This is an _independent_ dockerized version of the SIROCCO tools for easy deployment and sharing. A few of the dependencies are not directly available hence included in this repository under their respective LICENSE.

The main reason behind this docker is the SIROCCO tools only compiles with the older version of the GCC (GCC7), and does not play very well will new linux installations. Additionally, it does not allow using the tools in Windows or Mac. With docker, all these problems can be mitigated with extremely low overhead from docker itself (typically ~1MB), and with full access to the total computer capacity (processor, RAM) unlike virtualbox, or other virtualization systems. In other words, the SIROCCO tools can be used seamlessly in Windows, Linux, or Mac.

# Installation
To run the docker, one needs to install docker in their system. The docker is tested run on WINDOWS, LINUX, and MAC (both Intel and Apple Silicon chips). 

Step 1. Install Docker for [Windows](https://docs.docker.com/desktop/install/windows-install/), [Linux](https://docs.docker.com/desktop/install/linux-install/), [Mac](https://docs.docker.com/desktop/install/mac-install/). The installer size of docker is ~600MB.

Step 2. Open your terminal (linux, mac) or cmd/powershell (windows), and run `docker pull jamal919/tools`. This will download the docker which contains the compiled version of SIROCCO tools.

# Usage
The docker images are run using `docker run <imagename>` command, which creates a container from the image, run the tasks, and exit. However, for practical purposes, one needs to add some extra flags to the run command. 

## Interactive mode
To run the docker with mounting a certain location to `/mnt` of the docker in interactive mode (`-it`) and without retaining the exited docker (`--rm`) the following command can be used - 

`docker run -it --rm -v /path/in/the/local/system:/mnt jamal919/tools`

If you want to run the docker in the current directory, use the following command - 

```
# In linux or mac
docker run -it --rm -v `pwd`:/mnt jamal919/tools

# In windows
docker run -it --rm -v ${PWD}:/mnt jamal919/tools
```

Once inside the interactive terminal, one can access all the commands available from SIROCCO tools.

## Non-interactive mode
If you want to run the commands in non-interactive mode, e.g., the commands will be executed and the docker will exit, you can remove `-it` flag and provide your command at the end of the previously shown commands. For example, if you want to run `comodo-detidor` on a file located in the current directory `pwd`, then the command will look like the following - 

```
docker run -it --rm -v `pwd`:/mnt jamal919/tools comodo-detidor <option1> <option2> <file> <etc>
```