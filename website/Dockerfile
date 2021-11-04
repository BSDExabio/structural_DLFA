FROM python:3.8-alpine

RUN apt-get update -yqq && apt-get upgrade -yqq ca-certificates

RUN mkdir /app

WORKDIR /app

ADD requirements.txt /app
ADD main.py /app

RUN pip install --no-cache-dir --upgrade pip && \
      pip install --no-cache-dir -r requirements.txt

CMD ["gunicorn", "-w 4", "-b", "0.0.0.0:8000", "main:app"]
