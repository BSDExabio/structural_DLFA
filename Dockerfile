FROM python:3.8-slim

ENV FLASK_APP=/usr/src/website/app.py
ENV PYTHONPATH=/usr/src/

COPY website/requirements.txt ./
RUN pip install --no-cache-dir --upgrade pip && \
      pip install --prefer-binary --no-cache-dir -r requirements.txt

COPY website .
COPY database .

WORKDIR /usr/src/website
CMD ["flask", "run", "--host=0.0.0.0"]
