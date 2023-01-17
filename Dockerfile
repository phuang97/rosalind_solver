FROM python:3.9

WORKDIR /rosalind_solver
COPY requirements.txt ./
RUN apt update
RUN pip install --upgrade pip
RUN pip install --no-cache-dir \
    -r requirements.txt 
#COPY setup.py ./
#COPY rosalind_solver/__init__.py ./rosalind_solver/
COPY . ./
RUN python setup.py develop
