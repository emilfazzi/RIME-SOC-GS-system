FROM python:3.10


WORKDIR /usr/src

RUN apt-get update
RUN apt install -y python3-pip
RUN pip install --upgrade pip

COPY * ./

RUN pip3 install -r requirements.txt


RUN pip install --no-cache-dir -r requirements.txt



CMD ["python3", "tm2raw.py"]

