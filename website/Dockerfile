FROM python:3.8-alpine

RUN apk update && apk upgrade

RUN mkdir /app

WORKDIR /app

ADD requirements.txt .
ADD app.py .

RUN pip install --no-cache-dir --upgrade pip && \
      pip install --no-cache-dir -r requirements.txt

CMD ["gunicorn", "-w 4", "-b", "0.0.0.0:8000", "main:app"]
