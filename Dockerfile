FROM python:3.9


WORKDIR /rosalind_solver
COPY requirements.txt ./
COPY . ./
RUN apt update
RUN pip install --upgrade pip
RUN pip install --no-cache-dir \
    -r requirements.txt 
RUN python -m pip install -e .

