FROM python:3.9-slim

LABEL maintainer="colettima@ornl.gov"

ENV FLASK_APP=/usr/src/website/app.py
# FLASK_DEPLOYMENT can also be 'production'
ENV FLASK_DEPLOYMENT=development

ENV PYTHONPATH=/usr/src/

COPY website/requirements.txt ./
RUN pip install --no-cache-dir --upgrade pip && \
      pip install --prefer-binary --no-cache-dir -r requirements.txt

VOLUME dlfa_db

COPY ./db/dlfa.db dlfa_db/

COPY website/ /usr/src/website/
COPY database/ /usr/src/database/

WORKDIR /usr/src/website

CMD ["python3", "-m", "flask", "run"]
