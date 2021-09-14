# IDEPA Dockerfile
# Version 1.0

# Base ubuntu:16.04
# IDEPA Dockerfile
# Version 1.0


FROM continuumio/anaconda3 
WORKDIR /IDEPA-XMBD

# 修改源为国内源
RUN sed -i 's/archive.ubuntu.com/mirrors.aliyun.com/g' /etc/apt/sources.list
RUN sed -i 's/security.ubuntu.com/mirrors.aliyun.com/g' /etc/apt/sources.list

RUN apt-get clean
RUN apt-get update --fix-missing

# 更新系统
# 安装apt-utils（允许安装第三方软件）、sudo、vim、wget（下载文件）
RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get install -y apt-utils sudo vim wget openssh-server

RUN apt-get update
RUN apt-get install -y --no-install-recommends build-essential cmake git python-dev python-pip

# 修改ssh配置文件，文件位置/etc/ssh/sshd_config，添加允许所有主机连接，
RUN sed -i 's/#PermitRootLogin prohibit-password/PermitRootLogin yes/g' /etc/ssh/sshd_config
RUN echo 'sshd:ALL' >> /etc/hosts.aldlow
 
# Create the environment:
COPY . /IDEPA-XMBD/

RUN conda env create -f environment.yml
 
SHELL ["/bin/bash", "-c"]

 
RUN echo "source activate reoa3" > ~/.bashrc
RUN bash /IDEPA-XMBD/source_reoa3.sh
ENV PATH /opt/conda/envs/reoa3/bin:$PATH

RUN ["/opt/conda/bin/pip" ,"install", "rpy2==3.4.5"]
RUN python /IDEPA-XMBD/r_package_install.py

 
CMD ["/bin/bash"]

  
