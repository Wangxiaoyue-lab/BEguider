# Python 3.8-slim
FROM python:3.8-slim

#  wget  bzip2
RUN apt-get update && \
    apt-get install -y wget bzip2

#  Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda3-latest-Linux-x86_64.sh

# Setting environment variables

ENV PATH="/opt/conda/bin:${PATH}"

# environment.yml 
COPY environment.yml .

# conda 
RUN conda env create -f environment.yml

# env
ENV CONDA_DEFAULT_ENV=beguider
ENV CONDA_PREFIX=/opt/conda/envs/$CONDA_DEFAULT_ENV
ENV PATH=$CONDA_PREFIX/bin:$PATH

#  requirements.txt 
COPY requirements.txt .
RUN pip install -r requirements.txt

# port
EXPOSE 8501

# Copy all files in the current directory to the /app directory in the container
COPY .  app/
#COPY . .

# WORKDIR
WORKDIR /app

HEALTHCHECK CMD curl --fail http://localhost:8501/_stcore/health

# RUN
CMD ["sh", "-c", "cd /app/webs_main && streamlit run Intro.py --server.port=8501 --server.address=0.0.0.0"]
