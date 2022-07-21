FROM python:3.9-slim

LABEL maintainer="colettima@ornl.gov"

ENV FLASK_APP=/usr/src/website/app.py
ENV PYTHONPATH=/usr/src/

COPY website/requirements.txt ./
RUN pip install --no-cache-dir --upgrade pip && \
      pip install --prefer-binary --no-cache-dir -r requirements.txt

VOLUME dlfa_db

COPY ./db/dlfa.db dlfa_db/

COPY website/ /usr/src/website/
COPY database/ /usr/src/database/

WORKDIR /usr/src/website

CMD ["python3", "-m", "flask", "run", "--port=80", "--host=128.219.186.120"]

